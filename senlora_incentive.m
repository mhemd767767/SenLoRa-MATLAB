function schedule = senlora_incentive(activity, opts)
%SENLORA_INCENTIVE Phase 5: offline simulated adaptive packet requests.
%   An entropy-triggered request needs BOTH two consecutive walking windows
%   AND a FINITE, VALID measured sample entropy kappa > gamma. Undefined
%   entropy (NaN or Inf) NEVER counts as a measured threshold crossing.
%   Previous versions requested packets when ~isfinite(kappa), causing
%   apparent "entropy" requests at packets 9,15,21,27 before the configured
%   30-packet / 12-second minimum. This version fixes that bug.
%
%   This is NOT a closed-loop LoRaWAN implementation: request times are
%   proposals only and the extra packets are not simulated or decoded.
%   Gamma, motion classification, and cooldown require calibration. The
%   paper's gamma=0.01496 is scenario-specific, NOT a universal threshold.
%
%   Diagnostics: kappa_history, kappa_valid, kappa_reason, kappa_input_packets,
%   kappa_window_start_idx, event_indices, event_reason, request_at.
if nargin<2 || isempty(opts), opts=struct(); end
if ~isfield(opts,'gamma'), opts.gamma=[]; end
if ~isfield(opts,'burst_packets'), opts.burst_packets=5; end
if ~isfield(opts,'motion_window_packets'),opts.motion_window_packets=8;end
if ~isfield(opts,'consecutive_windows'),opts.consecutive_windows=2;end
if ~isfield(opts,'entropy_min_packets'),opts.entropy_min_packets=30;end
if ~isfield(opts,'entropy_min_duration_s'),opts.entropy_min_duration_s=12;end
if ~isfield(opts,'cooldown_packets'),opts.cooldown_packets=5;end
if ~isfield(opts,'max_bursts'),opts.max_bursts=4;end
if ~isfield(opts,'entropy_opts')
    opts.entropy_opts=struct('m',2,'r_fraction',0.2);
end
if ~isfield(opts,'max_gap_factor'),opts.max_gap_factor=2.5;end
if ~isfield(opts,'walking_threshold')
    opts.walking_threshold=activity.walking_threshold;
end
validateattributes(opts.burst_packets,{'numeric'},{'scalar','integer','positive'});
validateattributes(opts.motion_window_packets,{'numeric'}, ...
    {'scalar','integer','>=',3});
validateattributes(opts.consecutive_windows,{'numeric'}, ...
    {'scalar','integer','positive'});
validateattributes(opts.entropy_min_packets,{'numeric'}, ...
    {'scalar','integer','>=',5});
validateattributes(opts.entropy_min_duration_s,{'numeric'}, ...
    {'scalar','real','finite','nonnegative'});
validateattributes(opts.max_gap_factor,{'numeric'}, ...
    {'scalar','real','finite','>',1});
validateattributes(opts.cooldown_packets,{'numeric'}, ...
    {'scalar','integer','nonnegative'});
validateattributes(opts.max_bursts,{'numeric'}, ...
    {'scalar','integer','nonnegative'});
if ~isempty(opts.gamma)
    validateattributes(opts.gamma,{'numeric'}, ...
        {'scalar','real','nonnegative','finite'});
end
t=double(activity.times_s(:)); y=double(activity.amplitude(:));
n=numel(y);
if numel(t)~=n || any(~isfinite(t)) || any(~isfinite(y)) || ...
        (n>1 && any(diff(t)<=0))
    error('Activity must have equal-length finite amplitude and strictly increasing times_s.');
end
schedule = struct('times_s',t,'kappa_history',nan(n,1), ...
    'kappa_valid',false(n,1), ...
    'kappa_reason',{repmat({'not_enough_samples'},n,1)}, ...
    'kappa_input_packets',zeros(n,1), ...
    'kappa_window_start_idx',zeros(n,1), ...
    'walking_window',false(n,1),'walking_runs',zeros(n,1), ...
    'request_at',false(n,1),'requested_packet_count',zeros(n,1), ...
    'event_reason',{repmat({''},n,1)}, ...
    'total_extra_packets',0,'event_indices',zeros(0,1), ...
    'gamma',opts.gamma,'paper_gamma_reference_only',0.01496, ...
    'reason','', 'simulated_only',true);
if isempty(opts.gamma)
    schedule.reason='gamma_not_calibrated';
    return;
end
if ~isscalar(opts.walking_threshold) || ...
        ~isfinite(opts.walking_threshold) || opts.walking_threshold<=0
    schedule.reason='walking_threshold_not_calibrated';
    return;
end
run_count=0; bursts=0; last_event=-Inf;
for i=1:n
    if i>=opts.motion_window_packets
        w=y(i-opts.motion_window_packets+1:i);
        w_score=std(w)/max(abs(median(w)),eps);
        iswalking=w_score>opts.walking_threshold;
        schedule.walking_window(i)=iswalking;
        if iswalking, run_count=run_count+1; else, run_count=0; end
    end
    schedule.walking_runs(i)=run_count;
    if i>=opts.entropy_min_packets && ...
            t(i)-t(1)>=opts.entropy_min_duration_s
        % Never interpolate across a substantial missing-packet gap.
        % Use the most recent gap-free suffix, then enforce minimum
        % observations and duration ON THAT SUFFIX (not the full history).
        delta=diff(t(1:i)); typical_dt=median(delta);
        if isfinite(typical_dt) && typical_dt>0
            last_gap=find(delta>opts.max_gap_factor*typical_dt,1,'last');
            start_idx=1;
            if ~isempty(last_gap), start_idx=last_gap+1; end
            valid_count=i-start_idx+1;
            schedule.kappa_input_packets(i)=valid_count;
            schedule.kappa_window_start_idx(i)=start_idx;
            if valid_count>=opts.entropy_min_packets && ...
                    t(i)-t(start_idx)>=opts.entropy_min_duration_s
                sub_t=t(start_idx:i);
                dt=median(diff(sub_t));
                grid=(sub_t(1):dt:sub_t(end))';
                vals=interp1(sub_t,y(start_idx:i),grid,'linear');
                [value,detail]=senlora_sample_entropy(vals,opts.entropy_opts);
                schedule.kappa_history(i)=value;
                schedule.kappa_valid(i)=detail.valid && isfinite(value);
                schedule.kappa_reason{i}=detail.reason;
            else
                schedule.kappa_reason{i}='insufficient_contiguous_observations';
            end
        else
            schedule.kappa_reason{i}='invalid_sampling_interval';
        end
    end
    if run_count<opts.consecutive_windows || ...
            i-last_event<=opts.cooldown_packets || bursts>=opts.max_bursts
        continue;
    end
    % The critical fix: NaN/Inf/missing entropy cannot trigger an
    % entropy-based request. It is not evidence that kappa > gamma.
    if schedule.kappa_valid(i) && ...
            schedule.kappa_history(i)>opts.gamma
        schedule.request_at(i)=true;
        schedule.requested_packet_count(i)=opts.burst_packets;
        schedule.event_reason{i}='measured_kappa_above_gamma';
        last_event=i;
        bursts=bursts+1;
    end
end
schedule.event_indices=find(schedule.request_at);
schedule.total_extra_packets=sum(schedule.requested_packet_count);
schedule.reason='evaluated';
end
