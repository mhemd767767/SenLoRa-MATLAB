% BREATHING_18BPM SenLoRa respiration sensing at 18 BPM.
%
%   This experiment sends actual simulated LoRa IQ through the current
%   iterative soft receiver, reconstructs confidence-symbol CSI, and
%   estimates the respiration frequency from packet-mean CSI amplitude.
%   The injected breathing rate and channel are retained for validation.
%
%   This is a standalone breathing-activity test using simulated IQ.
%   No hard receiver methods are required.
%   The MATLAB simulation's m=2 and r=0.2*std entropy parameters are
%   assumptions, not published SenLoRa settings. This test focuses on BPM.
%
%   Run from the project root:
%      run('examples/breathing_18bpm.m')
%   Core files must be on the MATLAB path (see README.md).

clearvars; clc;

% -------------------------------------------------------------------------
% LoRa PHY and ground-truth breathing parameters.
% -------------------------------------------------------------------------
phy = LoRaPHY(470e6, 7, 125e3, 250e3);
phy.fast_mode = true;
phy.cr = 4;
phy.crc = 1;

truth_bpm = 18.0;
truth_hz = truth_bpm / 60;                    % 0.30 Hz

% One-second periodic sampling resolves a 0.30 Hz breathing tone:
% the packet-rate Nyquist frequency is 0.50 Hz.
opts = struct();
opts.mode = 'iq';
opts.adaptive = false;                       % Isolate respiration sensing
opts.duration_s = 45;                       % 46 received packets
opts.interval_s = 1.0;
opts.walk_end_s = 0;                         % Breathing starts at t = 0
opts.breath_end_s = opts.duration_s + 1;     % Breathing for the entire last packet
opts.breath_hz = truth_hz;
opts.breath_m = 0.008;                       % 8-mm displacement assumption
opts.snr_db = 20;
opts.snr_reference = 'received';             % Actual received SNR
opts.cfo_hz = 30;
opts.cr = 4;
opts.max_iter = 5;
opts.seed = 1818;
opts.max_bursts = 0;
opts.activity_min_packets = 12;
opts.activity_min_duration_s = 15;
opts.activity_min_r2 = 0.25;
opts.frequency_hz = (0.10:0.0025:0.45)';
opts.entropy_opts = struct('m',2,'r_fraction',0.2);
opts.throw_receiver_errors = true;
opts.csi_opts = struct('guard_samples',2, 'cfo_correction',true, ...
    'remove_static',true, 'iter_opts',struct('pilot_fallback',true));

fprintf('SenLoRa respiration test: true rate = %.2f BPM (%.3f Hz)\n', ...
    truth_bpm, truth_hz);
fprintf('SF=%d, BW=%.0f kHz, received SNR=%.1f dB, CFO=%.1f Hz\n', ...
    phy.sf,phy.bw/1e3,opts.snr_db,opts.cfo_hz);

% -------------------------------------------------------------------------
% Real signal path: TX chirps -> channel + noise -> iterative soft RX
% -> confidence-symbol CSI -> multi-packet respiration estimator.
% -------------------------------------------------------------------------
received = phy.senlora_phase6(opts);
act = received.activity;

fprintf('\nIQ packets: TX=%d, payload correct=%d, usable CSI=%d\n', ...
    received.transmitted_packets,received.payload_correct_packets, ...
    received.usable_packets);
fprintf('Selected confidence chirps: %d (mean %.2f per usable packet)\n', ...
    received.confidence_chirps, received.confidence_chirps / ...
    max(received.usable_packets,1));
fprintf('Soft IQ: estimated BPM %.3f, R^2 %.3f, valid=%d\n', ...
    act.respiration_bpm,act.respiration_r2,act.respiration_valid);
if isfinite(act.respiration_bpm)
    fprintf('Soft IQ absolute BPM error = %.3f BPM\n', ...
        abs(act.respiration_bpm-truth_bpm));
end

% -------------------------------------------------------------------------
% Compare extracted packet CSI amplitude to the same-timestamp channel.
% The reference channel is used only for evaluation, not decoding.
% -------------------------------------------------------------------------
t = act.times_s(:);
y = act.amplitude(:);
records = received.records;
usable = logical([records.usable_csi]);
truth_t = [records(usable).time_s]';
truth_a = [records(usable).truth_mean_amplitude]';
if ~isempty(t) && isequal(size(t),size(truth_t)) && ...
        max(abs(t-truth_t))<1e-8
    truth_aligned = truth_a;
else
    truth_aligned = interp1(truth_t,truth_a,t,'linear',NaN);
end
valid_pair = isfinite(truth_aligned) & isfinite(y);
amp_mae = NaN;
amp_corr = NaN;
if any(valid_pair)
    amp_mae = mean(abs(y(valid_pair)-truth_aligned(valid_pair)));
end
if nnz(valid_pair)>=3 && std(y(valid_pair))>0 && ...
        std(truth_aligned(valid_pair))>0
    rho = corrcoef(y(valid_pair),truth_aligned(valid_pair));
    amp_corr = rho(1,2);
end
fprintf('Packet CSI amplitude MAE=%.5f; correlation=%.4f\n', ...
    amp_mae,amp_corr);

% -------------------------------------------------------------------------
% Plot the observed CSI, injected channel truth, and frequency-search R^2.
% -------------------------------------------------------------------------
figure('Name','SenLoRa: 18 BPM breathing via soft LoRa CSI');
subplot(2,1,1);
plot(t,truth_aligned,'k--','LineWidth',1.4, ...
    'DisplayName','Injected true channel'); hold on;
plot(t,y,'b.','MarkerSize',12,'DisplayName','Soft receiver packet CSI');
if act.respiration_valid
    plot(t,act.respiration_fit,'r-','LineWidth',1.5, ...
        'DisplayName','Best sinusoidal fit to soft CSI');
end
grid on;
xlabel('Time (s)'); ylabel('Mean |CSI| / packet');
title(sprintf('Breathing truth %.1f BPM; estimated %.2f BPM; R^2=%.3f', ...
    truth_bpm,act.respiration_bpm,act.respiration_r2));
legend('Location','best');

subplot(2,1,2);
freq = opts.frequency_hz(:);
score = local_frequency_r2(t,y,freq);
plot(60*freq,score,'b-','LineWidth',1.5, ...
    'DisplayName','Soft CSI sinusoid R^2'); hold on;
xline(truth_bpm,'k--','True 18 BPM','LineWidth',1.3);
if isfinite(act.respiration_bpm)
    xline(act.respiration_bpm,'r:','Estimated BPM','LineWidth',1.3);
end
grid on;
xlabel('Candidate breathing frequency (BPM)');
ylabel('Sine-fit R^2');
title('Respiration frequency search (packet timestamps)');
legend('Location','best');

% -------------------------------------------------------------------------
% Validate the actual soft IQ path against injected 18 BPM truth.
% Never substitute the injected rate when receiver sensing fails.
% -------------------------------------------------------------------------
assert(received.payload_correct_packets>0 && received.usable_packets>=12, ...
    'Too few correctly received packets / usable CSI for respiration.');
if act.respiration_valid && abs(act.respiration_bpm-truth_bpm)<=1.0
    fprintf('PASS: soft IQ receiver recovered 18 BPM within 1 BPM.\n');
else
    warning('18 BPM NOT validated on soft IQ: %.3f BPM (R^2 %.3f).', ...
        act.respiration_bpm,act.respiration_r2);
end

function score=local_frequency_r2(t,y,freq)
%LOCAL_FREQUENCY_R2 Scan sinusoidal fits without rank-deficient candidates.
t=t(:); y=y(:); freq=freq(:);
score=nan(size(freq));
if numel(y)<3,return;end
ss0=sum((y-mean(y)).^2);
if ~isfinite(ss0) || ss0<=eps,return;end
for k=1:numel(freq)
    X=[ones(numel(t),1),cos(2*pi*freq(k)*(t-t(1))), ...
        sin(2*pi*freq(k)*(t-t(1)))];
    if rank(X)<3,continue;end
    fit=X*(X\y);
    score(k)=1-sum((y-fit).^2)/ss0;
end
end
