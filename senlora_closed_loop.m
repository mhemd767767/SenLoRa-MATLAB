function result = senlora_closed_loop(phy, opts)
%SENLORA_CLOSED_LOOP Phase 6: event-driven, closed-loop SenLoRa simulation.
%
%   RESULT = SENLORA_CLOSED_LOOP(PHY, OPTS) generates periodic LoRa
%   transmissions, demodulates their IQ using the iterative soft receiver,
%   extracts confidence-symbol CSI, and injects ACTUAL extra simulated
%   transmissions when the simulated controller requests them.
%
%   Unlike the Phase 5 scheduler, every requested packet is synthesized,
%   demodulated and incorporated into subsequent entropy measurements.
%   The packet-event queue is updated only AFTER observing each packet.
%
%   Paper parameters: two consecutive walking detections, five consecutive
%   incentive packets, gamma=0.01496 (SenLoRa Sec. 7, Fig. 8a).
%   The paper's threshold is calibrated on a different CSI representation;
%   using that number here is a REFERENCE test, not a reproduction of its
%   calibrated respiration decision boundary.
%
%   Simulation assumptions, NOT specified by the paper: packet interval,
%   packet burst spacing, m=2, r=0.2*std, 0.035 walking threshold,
%   four-packet motion window, AWGN, one-path moving reflector, cooldown,
%   and max burst limit. See README.md for parameter assumptions.
%
%   opts.mode = 'iq' (only mode); opts.adaptive = true/false.
%   opts.walk_end_s and opts.breath_end_s define motion transitions.
%
%   The separate bootstrap trigger after walking is intentionally labeled
%   'walking_bootstrap': it does NOT claim to measure unavailable entropy.
%   Later requests require finite measured kappa > gamma.
if nargin<2 || isempty(opts), opts=struct(); end
opts = local_defaults(opts);
state = rng;
cleanup_rng = onCleanup(@() rng(state)); %#ok<NASGU>
validateattributes(opts.gamma,{'numeric'}, ...
    {'scalar','finite','nonnegative'});
validateattributes(opts.burst_packets,{'numeric'}, ...
    {'scalar','integer','positive'});
validateattributes(opts.interval_s,{'numeric'}, ...
    {'scalar','positive','finite'});
validateattributes(opts.burst_spacing_s,{'numeric'}, ...
    {'scalar','positive','finite'});
if ~strcmp(opts.mode,'iq')
    error('mode must be iq.');
end
if ~any(strcmp(opts.snr_reference,{'transmitted','received'}))
    error('snr_reference must be transmitted or received.');
end
if ~isscalar(opts.breath_hz) || opts.breath_hz<=0
    error('breath_hz must be a positive scalar.');
end
if ~isscalar(opts.duration_s) || opts.duration_s<=0
    error('duration_s must be positive.');
end
if opts.walk_end_s >= opts.breath_end_s
    error('walk_end_s must precede breath_end_s.');
end
if opts.burst_spacing_s <= opts.packet_guard_s
    error('burst spacing must exceed packet guard.');
end

% Build one repeatable transmitted payload and chirp waveform.
phyTx = LoRaPHY(phy.rf_freq,phy.sf,phy.bw,phy.fs);
phyTx.cr=opts.cr;
phyTx.crc=1;
phyTx.fast_mode=true;
payload=uint8(opts.payload(:));
symbols=phyTx.encode(payload);
tx=phyTx.modulate(symbols);
packet_duration_s=numel(tx)/phy.fs;
min_separation = packet_duration_s + opts.packet_guard_s;
if opts.burst_spacing_s<min_separation
    error('Burst spacing %.3fs must exceed packet airtime plus guard %.3fs.', ...
        opts.burst_spacing_s,min_separation);
end
if opts.interval_s<min_separation
    error('Baseline packets overlap; increase opts.interval_s.');
end

% Store timestamps and their source. 0=periodic, 1=incentive.
event_time=(0:opts.interval_s:opts.duration_s)';
event_type=zeros(numel(event_time),1);
base_count=numel(event_time);
record=struct('time_s',{},'is_extra',{},'usable_csi',{}, ...
    'payload_correct',{},'pilot_count',{},'truth_mean_amplitude',{},'kappa',{}, ...
    'kappa_valid',{},'kappa_reason',{},'motion_score',{},'walking',{}, ...
    'event_reason',{},'added_packets',{});
csi=struct([]);
chirp_amplitude=zeros(0,1);
chirp_time=zeros(0,1);
truth_chirp_amplitude=zeros(0,1);
request_times=zeros(0,1);
request_reasons=cell(0,1);
request_kappa=zeros(0,1);
request_numbers=zeros(0,1);
bootstrap_done=false;
session_active=false;
previous_walk=false;
walk_run=0;
last_request=-Inf;

% Process packets in chronological order: cannot inspect future packets.
i=1;
while i<=numel(event_time)
    t=event_time(i);
    is_extra=logical(event_type(i));
    [packet,decoded_ok] = local_receive(phy,opts,tx,payload,t);
    usable=~isempty(packet) && ~isempty(packet.corrected_csi) && ...
        all(isfinite(real(packet.corrected_csi(:)))) && ...
        all(isfinite(imag(packet.corrected_csi(:))));
    pilots=0;
    truth_mean_amplitude=NaN;
    if usable
        packet.packet_time_s=t;
        if isfield(packet,'pilot_times_s') && ~isempty(packet.pilot_times_s)
            pt=t+packet.pilot_times_s(:);
        else
            pt=t+(1:numel(packet.corrected_csi))'*0.001;
        end
        h=packet.corrected_csi(:);
        chirp_amplitude=[chirp_amplitude;abs(h)]; %#ok<AGROW>
        chirp_time=[chirp_time;pt]; %#ok<AGROW>
        % Ground truth is generated ONLY for simulation diagnostics.
        % This is not a channel estimate available to a real gateway.
        true_amp=abs(local_channel(pt,opts,phy.rf_freq));
        truth_chirp_amplitude=[truth_chirp_amplitude;true_amp]; %#ok<AGROW>
        truth_mean_amplitude=mean(true_amp);
        pilots=numel(h);
        if isempty(csi)
            csi=packet;
        else
            csi(end+1)=packet; %#ok<AGROW>
        end
    end
    nCSI=numel(csi);

    % Chirp-index sample entropy as in Fig. 4's x-axis.  The physical
    % packet gaps are not suppressed; report this interpretation explicitly.
    % This is not interchangeable with packet-mean sample entropy.
    kappa=NaN; kappa_valid=false;
    kappa_reason='not_enough_confidence_chirps';
    if numel(chirp_amplitude)>=opts.min_entropy_chirps
        kappa_reason='awaiting_usable_csi_or_update';
    end
    if numel(chirp_amplitude)>=opts.min_entropy_chirps && ...
            mod(nCSI,opts.entropy_update_packets)==0 && usable
        idx0=max(1,numel(chirp_amplitude)-opts.entropy_max_chirps+1);
        a=chirp_amplitude(idx0:end);
        [kappa,detail]=senlora_sample_entropy(a,opts.entropy_opts);
        kappa_valid=detail.valid && isfinite(kappa);
        kappa_reason=detail.reason;
    end

    motion=NaN; walking=false;
    if nCSI>=opts.motion_window_packets
        idx0=max(1,nCSI-opts.motion_window_packets+1);
        recent=csi(idx0:end);
        amps=zeros(numel(recent),1);
        for j=1:numel(recent)
            amps(j)=mean(abs(recent(j).corrected_csi(:)));
        end
        motion=std(amps)/max(abs(median(amps)),eps);
        walking=motion>opts.walking_threshold;
    end
    if usable && walking
        if previous_walk, walk_run=walk_run+1; else, walk_run=1; end
    elseif usable
        walk_run=0;
    end
    previous_walk=walking;

    reason=''; nAdded=0;
    if opts.adaptive && usable && ...
            t-last_request >= opts.cooldown_s && ...
            numel(request_times)<opts.max_bursts
        if ~bootstrap_done && walk_run>=opts.consecutive_windows
            % SenLoRa Sec. 7: activity-triggered bootstrap after two
            % positive observations. No entropy assertion is made here.
            reason='walking_bootstrap';
        elseif session_active && kappa_valid && kappa>opts.gamma
            % Stay active during the sensing session even after walking
            % stops: otherwise additional breathing packets cannot fire.
            reason='finite_entropy_above_gamma';
        end
    end
    if ~isempty(reason)
        proposed=t+(1:opts.burst_packets)'*opts.burst_spacing_s;
        if proposed(end)<=opts.duration_s
            % Avoid overlaps with existing scheduled transmissions.
            for j=1:numel(proposed)
                while any(abs([event_time;proposed(1:j-1)]-proposed(j)) ...
                        <min_separation)
                    proposed(j)=proposed(j)+min_separation;
                end
            end
            if proposed(end)<=opts.duration_s
                event_time=[event_time;proposed]; %#ok<AGROW>
                event_type=[event_type;ones(numel(proposed),1)]; %#ok<AGROW>
                [event_time,ord]=sort(event_time);
                event_type=event_type(ord);
                nAdded=numel(proposed);
                request_times(end+1,1)=t; %#ok<AGROW>
                request_reasons{end+1,1}=reason; %#ok<AGROW>
                request_kappa(end+1,1)=kappa; %#ok<AGROW>
                request_numbers(end+1,1)=nAdded; %#ok<AGROW>
                last_request=t;
                bootstrap_done=true;
                session_active=true;
            else
                reason='';
            end
        else
            reason='';
        end
    end
    record(end+1)=struct('time_s',t,'is_extra',is_extra, ...
        'usable_csi',usable,'payload_correct',decoded_ok, ...
        'pilot_count',pilots, ...
        'truth_mean_amplitude',truth_mean_amplitude, ...
        'kappa',kappa,'kappa_valid',kappa_valid, ...
        'kappa_reason',kappa_reason, ...
        'motion_score',motion,'walking',walking, ...
        'event_reason',reason,'added_packets',nAdded); %#ok<AGROW>
    i=i+1;
end

result=struct();
result.records=record;
result.csi=csi;
result.chirp_amplitude=chirp_amplitude;
result.chirp_time_s=chirp_time;
result.truth_chirp_amplitude=truth_chirp_amplitude;
result.request_times_s=request_times;
result.request_reasons=request_reasons;
result.request_kappa=request_kappa;
result.request_numbers=request_numbers;
result.sent_times_s=[record.time_s]';
result.sent_is_extra=logical([record.is_extra]');
result.base_packets=base_count;
result.extra_packets=sum(result.sent_is_extra);
result.transmitted_packets=numel(record);
result.usable_packets=sum([record.usable_csi]);
result.payload_correct_packets=sum([record.payload_correct]);
% A decoded payload does not imply at least one confidence-symbol pilot.
result.zero_pilot_packets=sum([record.pilot_count]==0);
result.confidence_chirps=sum([record.pilot_count]);
result.valid_entropy_points=sum([record.kappa_valid]);
result.gamma=opts.gamma;
result.mode=opts.mode;
result.options=opts;
result.paper_parameters=struct('gamma',0.01496, ...
    'burst_packets',5,'walking_consecutive',2,'alpha',0.8);
result.is_hardware_replication=false;
% Consecutive packet CSI at actual event times: do not fabricate missing CSI.
activityopts=struct('min_packets',opts.activity_min_packets, ...
    'min_duration_s',opts.activity_min_duration_s, ...
    'min_fit_r2',opts.activity_min_r2, ...
    'frequency_hz',opts.frequency_hz);
if isempty(csi)
    result.activity=struct('times_s',[],'amplitude',[], ...
        'respiration_bpm',NaN,'respiration_r2',NaN,'respiration_valid',false);
else
    % Limit the fit to the breathing segment, not the walking episode.
    breathCSI=csi([csi.packet_time_s]>=opts.walk_end_s & ...
        [csi.packet_time_s]<=opts.breath_end_s);
    if isempty(breathCSI), breathCSI=csi; end
    result.activity=senlora_activity(breathCSI,activityopts);
end
end

function [packet,payload_ok] = local_receive(phy,opts,tx,payload,t)
%LOCAL_RECEIVE Generate one IQ packet and actually pass through soft RX.
packet=[]; payload_ok=false;
rng(opts.seed+round(t*1000),'twister');
n=numel(tx);
sample_times=t+(0:n-1)'/phy.fs;
h=local_channel(sample_times,opts,phy.rf_freq);
% The original Phase 6 code defined SNR relative to unit TX amplitude.
% Optionally define it relative to the actual received signal power.
% RX normalization is a SIMULATION CHOICE, not an author parameter.
noise_power=10^(-opts.snr_db/10);
if strcmp(opts.snr_reference,'received')
    noise_power=noise_power*mean(abs(h.*tx).^2);
end
noise_std=sqrt(noise_power/2);
w=noise_std*(randn(n,1)+1i*randn(n,1));
rx=h.*tx.*exp(1i*2*pi*opts.cfo_hz*sample_times)+w;
receiver=LoRaPHY(phy.rf_freq,phy.sf,phy.bw,phy.fs);
receiver.fast_mode=true;
receiver.cr=opts.cr;
receiver.crc=1;
try
    [out,decoded]=receiver.extract_csi(rx,opts.max_iter,opts.csi_opts);
    if ~isempty(out), packet=out(1); end
    if ~isempty(decoded) && size(decoded,1)>=numel(payload)
        payload_ok=isequal(uint8(decoded(1:numel(payload),1)),payload);
    end
catch ex
    if opts.throw_receiver_errors
        rethrow(ex);
    end
end
end

function h=local_channel(t,opts,rf)
%LOCAL_CHANNEL Static path plus one moving reflected propagation path.
lambda=299792458/rf;
d=opts.distance_m*ones(size(t));
walking=t<opts.walk_end_s;
breathing=t>=opts.walk_end_s & t<=opts.breath_end_s;
d(walking)=d(walking)+opts.walk_m*sin(2*pi*opts.walk_hz*t(walking));
d(breathing)=d(breathing)+opts.breath_m* ...
    sin(2*pi*opts.breath_hz*(t(breathing)-opts.walk_end_s));
h=opts.static_path+opts.reflected_path*exp(-1i*4*pi*d/lambda);
end

function opts=local_defaults(opts)
%LOCAL_DEFAULTS Paper-grounded parameters and explicit simulation choices.
def=struct('mode','iq','adaptive',true,'seed',2026, ...
    'duration_s',72,'interval_s',3,'burst_spacing_s',0.50, ...
    'packet_guard_s',0.02,'cr',4,'payload',uint8((1:12)'), ...
    'max_iter',3,'csi_opts',struct('guard_samples',2, ...
        'cfo_correction',true,'remove_static',true, ...
        'iter_opts',struct('pilot_fallback',true)), ...
    'snr_db',10,'snr_reference','transmitted', ...
    'cfo_hz',30,'breath_hz',0.25, ...
    'walk_hz',0.40,'breath_m',0.008,'walk_m',0.030, ...
    'walk_end_s',14,'breath_end_s',69, ...
    'distance_m',0.14,'static_path',0.85+0.05i, ...
    'reflected_path',0.45, ...
    'gamma',0.01496,'burst_packets',5, ...
    'consecutive_windows',2,'motion_window_packets',4, ...
    'walking_threshold',0.035,'cooldown_s',5,'max_bursts',4, ...
    'min_entropy_chirps',25,'entropy_max_chirps',600, ...
    'entropy_update_packets',1, ...
    'entropy_opts',struct('m',2,'r_fraction',0.2), ...
    'activity_min_packets',8,'activity_min_duration_s',12, ...
    'activity_min_r2',0.25, ...
    'frequency_hz',(0.10:0.0025:0.65)', ...
    'throw_receiver_errors',true);
f=fieldnames(def);
for j=1:numel(f)
    if ~isfield(opts,f{j}),opts.(f{j})=def.(f{j});end
end
if opts.max_bursts<0 || opts.cooldown_s<0 || ...
        opts.motion_window_packets<2 || opts.entropy_update_packets<1
    error('Invalid triggering parameters.');
end
end
