function out = senlora_activity(csi, opts)
%SENLORA_ACTIVITY Phase 4: multi-packet CSI motion and respiration sensing.
%   OUT = SENLORA_ACTIVITY(CSI, OPTS) uses complex CSI returned by the
%   soft-only LoRaPHY.extract_csi(). The feature is the mean *magnitude*
%   of corrected confidence-pilot CSI in each packet. Absolute phase is
%   deliberately NOT combined across packets: Phase 3 provides only
%   packet-relative CFO correction, not cross-packet phase calibration.
%
%   Respiration: scan sinusoidal frequencies via least-squares fit on
%   actual packet timestamps, robust to unequal inter-packet intervals.
%   Detection thresholds and sinusoidal fitting are SIMULATION HEURISTICS,
%   not specified classifiers in the published SenLoRa paper.
%
%   OPTS.stationary_csi: a stationary calibration capture (recommended).
%   OPTS.motion_threshold / walking_threshold: alternatively fixed scores.
%   OPTS.min_duration_s (12), min_cycles (2.5), min_packets (30).
%   OPTS.frequency_hz = 0.10:0.0025:0.65; OPTS.min_fit_r2=0.25.
if nargin < 2 || isempty(opts), opts = struct(); end
if ~isfield(opts,'min_duration_s'), opts.min_duration_s = 12; end
if ~isfield(opts,'min_cycles'), opts.min_cycles = 2.5; end
if ~isfield(opts,'min_packets'), opts.min_packets = 30; end
if ~isfield(opts,'frequency_hz')
    opts.frequency_hz = (0.10:0.0025:0.65)';
end
if ~isfield(opts,'min_fit_r2'), opts.min_fit_r2 = 0.25; end
if ~isfield(opts,'entropy_opts') || isempty(opts.entropy_opts)
    opts.entropy_opts = struct('m',2,'r_fraction',0.2);
end
if ~isfield(opts.entropy_opts,'m'), opts.entropy_opts.m=2; end
if ~isfield(opts,'max_gap_factor'), opts.max_gap_factor = 2.5; end
if ~isfield(opts,'motion_threshold'), opts.motion_threshold = NaN; end
if ~isfield(opts,'walking_threshold'), opts.walking_threshold = NaN; end
if ~isstruct(csi), error('CSI must be a struct array.'); end
[t,y,pilots,usable_ids] = local_csi_feature(csi);
num = numel(y);
span = 0;
if num >= 2, span = t(end)-t(1); end
out = struct('times_s',t,'amplitude',y, ...
    'pilot_counts',pilots,'usable_csi_indices',usable_ids, ...
    'usable_packets',num,'duration_s',span, ...
    'motion_score',NaN,'motion_threshold',NaN, ...
    'motion_detected',NaN,'walking_threshold',NaN, ...
    'walking_detected',NaN,'respiration_hz',NaN, ...
    'respiration_bpm',NaN,'respiration_r2',NaN, ...
    'respiration_valid',false,'respiration_fitted_amplitude',NaN, ...
    'respiration_fit',nan(num,1),'kappa',NaN, ...
    'entropy_valid',false,'entropy_detail',struct(), ...
    'entropy_sampling_note','Packet-mean CSI; approximate uniform-spacing required');
if num < 3, return; end
out.motion_score = std(y)/max(abs(median(y)),eps);
threshold = opts.motion_threshold;
walk_threshold = opts.walking_threshold;
if isfield(opts,'stationary_csi') && ~isempty(opts.stationary_csi)
    [~,quiet] = local_csi_feature(opts.stationary_csi);
    if numel(quiet)>=3
        baseline = std(quiet)/max(abs(median(quiet)),eps);
        if ~isfinite(threshold)
            threshold = max(4*baseline,0.002);
        end
        if ~isfinite(walk_threshold)
            walk_threshold = max(10*baseline,0.035);
        end
    end
end
out.motion_threshold = threshold;
out.walking_threshold = walk_threshold;
if isfinite(threshold), out.motion_detected = out.motion_score>threshold; end
if isfinite(walk_threshold)
    out.walking_detected = out.motion_score>walk_threshold;
end
if num < opts.min_packets || span < opts.min_duration_s, return; end
% Sample entropy is a sequence-index statistic, not naturally a
% nonuniform-time statistic. Regularize only if no substantial gap exists.
dt = median(diff(t));
if isfinite(dt) && dt>0 && ...
        all(diff(t) <= opts.max_gap_factor*dt)
    grid = (t(1):dt:t(end))';
    series = interp1(t,y,grid,'linear');
    if numel(series) >= opts.entropy_opts.m+3
        [out.kappa,out.entropy_detail] = ...
            senlora_sample_entropy(series,opts.entropy_opts);
        out.entropy_valid = out.entropy_detail.valid;
    end
else
    out.entropy_sampling_note = ...
        'Long packet gap: sample entropy not calculated';
end
% Fit a constant plus sine and cosine; actual times, no fabricated phase.
freqs = double(opts.frequency_hz(:));
if isempty(freqs) || any(~isfinite(freqs) | freqs<=0)
    error('frequency_hz must contain positive finite frequencies.');
end
sse0 = sum((y-mean(y)).^2);
if sse0 <= 1e-18, return; end
best_r2 = -Inf;
best_freq = NaN;
best_fit = nan(num,1);
best_amp = NaN;
for fi = 1:numel(freqs)
    f = freqs(fi);
    if f*span < opts.min_cycles, continue; end
    X = [ones(num,1),cos(2*pi*f*(t-t(1))),sin(2*pi*f*(t-t(1)))];
    % Regular packet times can alias some trial frequencies, leaving the
    % sine/cosine design matrix rank deficient (e.g., 0.25 Hz at 2 s).
    % These frequencies cannot identify all three model coefficients.
    % Do not fit them; nonuniform incentive times may restore full rank.
    if rank(X) < size(X,2), continue; end
    beta = X\y;
    fitted = X*beta;
    r2 = 1-sum((y-fitted).^2)/sse0;
    if r2>best_r2
        best_r2 = r2;
        best_freq = f;
        best_fit = fitted;
        best_amp = hypot(beta(2),beta(3));
    end
end
if ~isfinite(best_freq), return; end
out.respiration_hz = best_freq;
out.respiration_bpm = 60*best_freq;
out.respiration_r2 = best_r2;
out.respiration_fitted_amplitude = best_amp;
out.respiration_fit = best_fit;
out.respiration_valid = best_r2 >= opts.min_fit_r2;
end

function [times,amp,counts,source_ids] = local_csi_feature(csi)
times = zeros(0,1); amp = zeros(0,1); counts=zeros(0,1);
source_ids=zeros(0,1);
for k=1:numel(csi)
    item = csi(k);
    if ~isfield(item,'packet_time_s') || ...
            ~isfield(item,'corrected_csi')
        error('CSI entries must have packet_time_s and corrected_csi.');
    end
    v = item.corrected_csi(:);
    v = v(isfinite(real(v)) & isfinite(imag(v)));
    ts = double(item.packet_time_s);
    if isempty(v) || ~isscalar(ts) || ~isfinite(ts), continue; end
    times(end+1,1) = ts; %#ok<AGROW>
    amp(end+1,1) = mean(abs(v)); %#ok<AGROW>
    counts(end+1,1) = numel(v); %#ok<AGROW>
    source_ids(end+1,1) = k; %#ok<AGROW>
end
[times,ord] = sort(times);
amp = amp(ord); counts = counts(ord); source_ids=source_ids(ord);
if numel(times)>1 && any(diff(times)<=0)
    error('CSI packet timestamps must be strictly increasing.');
end
end
