% ADAPTIVE_PACKET_REQUESTS Closed-loop additional-packet simulation.
%
%   A simulated person walks, then breathes at 18 BPM. Compare the same
%   periodic LoRa packets with an adaptive controller. Additional packets
%   are generated as IQ, decoded with soft LoRa, and added to CSI only
%   when the simulated controller actually requests them.
%
%   This is a MATLAB-only demonstration: there are no LoRaWAN beacons.
%   Sample-entropy m=2, r=0.2*std, channel paths, and walking threshold
%   are simulation choices. The paper's experimental gamma=0.01496 is
%   retained for reference, but it is NOT calibrated for this model.

clear; clc; close all;
phy = LoRaPHY(470e6,7,125e3,250e3);
phy.fast_mode = true;
phy.cr = 4;
phy.crc = 1;

opts = struct('mode','iq','duration_s',45,'interval_s',1.0, ...
    'walk_end_s',12,'breath_end_s',46,'breath_hz',18/60, ...
    'snr_db',-4,'snr_reference','received','cfo_hz',30, ...
    'max_iter',5,'seed',2026,'gamma',0.01496, ...
    'burst_packets',5,'consecutive_windows',2, ...
    'burst_spacing_s',0.50,'max_bursts',3, ...
    'walking_threshold',0.025,'entropy_opts', ...
    struct('m',2,'r_fraction',0.2), ...
    'throw_receiver_errors',true);
fprintf('Adaptive packet requests: received SNR = %.1f dB\n', opts.snr_db);

% Baseline: periodic packet IQ processing, without incentive requests.
base_opts = opts;
base_opts.adaptive = false;
base = phy.senlora_phase6(base_opts);

% Adaptive: process the same baseline transmissions and any new IQ packets.
adaptive_opts = opts;
adaptive_opts.adaptive = true;
adaptive = phy.senlora_phase6(adaptive_opts);

fprintf('Periodic: TX=%d, decoded payloads=%d, usable CSI=%d\n', ...
    base.transmitted_packets,base.payload_correct_packets,base.usable_packets);
fprintf('Adaptive: TX=%d, decoded payloads=%d, usable CSI=%d\n', ...
    adaptive.transmitted_packets,adaptive.payload_correct_packets, ...
    adaptive.usable_packets);
fprintf('Extra packets actually processed: %d\n',adaptive.extra_packets);
fprintf('Requests issued: %d\n',numel(adaptive.request_times_s));
for k=1:numel(adaptive.request_times_s)
    fprintf('  %.2fs: %s, requested %d IQ packets, kappa %.4g\n', ...
        adaptive.request_times_s(k),adaptive.request_reasons{k}, ...
        adaptive.request_numbers(k),adaptive.request_kappa(k));
end
fprintf('18 BPM truth: periodic %.2f BPM, adaptive %.2f BPM\n', ...
    base.activity.respiration_bpm,adaptive.activity.respiration_bpm);
if adaptive.extra_packets==0
    fprintf(['No extra packets were triggered for this channel/threshold.\n', ...
        'This is a valid controller outcome, not evidence of adaptive gain.\n']);
end

figure('Name','Adaptive SenLoRa packet reception');
subplot(2,1,1);
plot(base.sent_times_s,zeros(size(base.sent_times_s)),'k.','MarkerSize',14);
hold on;
plot(adaptive.sent_times_s(adaptive.sent_is_extra), ...
    ones(nnz(adaptive.sent_is_extra),1),'r.','MarkerSize',16);
ylim([-0.5 1.5]); yticks([0 1]); yticklabels({'Periodic','Extra'});
xlabel('Packet time (s)');
title(sprintf('Packet schedule (received SNR %.1f dB)', opts.snr_db)); grid on;
subplot(2,1,2);
plot(base.activity.times_s,base.activity.amplitude,'o-', ...
    'DisplayName','Periodic CSI'); hold on;
plot(adaptive.activity.times_s,adaptive.activity.amplitude,'.-', ...
    'DisplayName','Adaptive CSI');
grid on; legend('Location','best');
xlabel('Packet time (s)'); ylabel('Mean |CSI|');
title(sprintf('Respiration sensing (received SNR %.1f dB)', opts.snr_db));

assert(adaptive.transmitted_packets==base.transmitted_packets+ ...
    adaptive.extra_packets,'Packet insertion count is inconsistent.');
assert(adaptive.extra_packets==sum(adaptive.request_numbers), ...
    'Requested packet accounting is inconsistent.');
