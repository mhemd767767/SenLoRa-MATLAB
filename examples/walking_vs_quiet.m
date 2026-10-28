% WALKING_VS_QUIET Compare CSI activity under static and moving channels.
%
%   This example sends two groups of LoRa IQ packets through the soft
%   receiver. The groups use the same noise seeds and packet timestamps,
%   but different physical channel motion. A heuristic amplitude-motion
%   score distinguishes walking from the quiet environment.
%
%   The walking threshold is an illustrative simulation choice; it is
%   NOT a published numerical parameter from the SenLoRa paper.

clear; clc; close all;
phy = LoRaPHY(470e6,7,125e3,250e3);
phy.fast_mode = true;
phy.cr = 4;
phy.crc = 1;

opts = struct('mode','iq','adaptive',false,'duration_s',30, ...
    'interval_s',1.0,'snr_db',-4,'snr_reference','received', ...
    'cfo_hz',30,'max_iter',5,'seed',2026,'max_bursts',0, ...
    'throw_receiver_errors',true);
fprintf('Walking vs quiet: received SNR = %.1f dB\n', opts.snr_db);

% The quiet channel remains static after the first observation.
quiet_opts = opts;
quiet_opts.walk_end_s = 0;
quiet_opts.breath_end_s = 0.001;
quiet = phy.senlora_phase6(quiet_opts);

% A 0.40-Hz reflector displacement is present throughout walking run.
walk_opts = opts;
walk_opts.walk_end_s = opts.duration_s+1;
walk_opts.breath_end_s = opts.duration_s+2;
walk = phy.senlora_phase6(walk_opts);

motion_opts = struct('walking_threshold',0.035,'motion_threshold',0.010, ...
    'min_packets',8,'min_duration_s',8);
quiet_act = phy.sense_activity(quiet.csi,motion_opts);
walk_act = phy.sense_activity(walk.csi,motion_opts);

fprintf('Quiet: IQ packets %d, correct payloads %d, usable CSI %d\n', ...
    quiet.transmitted_packets,quiet.payload_correct_packets,quiet.usable_packets);
fprintf('Walk:  IQ packets %d, correct payloads %d, usable CSI %d\n', ...
    walk.transmitted_packets,walk.payload_correct_packets,walk.usable_packets);
fprintf('Quiet motion score: %.4f; walking detected: %g\n', ...
    quiet_act.motion_score,quiet_act.walking_detected);
fprintf('Walk motion score:  %.4f; walking detected: %g\n', ...
    walk_act.motion_score,walk_act.walking_detected);

figure('Name','Walking versus quiet — decoded LoRa CSI');
plot(quiet_act.times_s,quiet_act.amplitude,'o-','LineWidth',1.2); hold on;
plot(walk_act.times_s,walk_act.amplitude,'o-','LineWidth',1.2);
grid on;
xlabel('Packet time (s)'); ylabel('Mean confidence-chirp |CSI|');
legend('Quiet','Walking','Location','best');
title(sprintf('Human activity from soft LoRa CSI (SNR %.1f dB)', opts.snr_db));

assert(quiet.usable_packets>=10 && walk.usable_packets>=10, ...
    'Insufficient usable CSI packets.');
if isfinite(walk_act.motion_score) && isfinite(quiet_act.motion_score)
    fprintf('Walking/quiet motion-score ratio: %.2f\n', ...
        walk_act.motion_score/max(quiet_act.motion_score,eps));
end
