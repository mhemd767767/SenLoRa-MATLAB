% COMMUNICATION_AND_CSI Decode one LoRa packet and reconstruct sensing CSI.
%
%   This example demonstrates the integrated communication and sensing
%   receive path using generated IQ and the iterative soft decoder.
%   The channel, CFO, noise and payload are controlled MATLAB inputs.
%   No hard LoRa decoding functions are used.

clear; clc; close all;
rng(2026,'twister');

% Configure the simulated LoRa transmitter and receiver.
rf = 470e6;
sf = 7;
bw = 125e3;
fs = 250e3;
phyTx = LoRaPHY(rf,sf,bw,fs);
phyTx.fast_mode = true;
phyTx.cr = 4;
phyTx.crc = 1;

payload = uint8((1:12)');
tx = phyTx.modulate(phyTx.encode(payload));

% Transmit through a single static path, oscillator offset, and AWGN.
snr_db = -4;
cfo_hz = 30;
h = 0.85 + 0.10i;
fprintf('Communication + CSI test: received SNR = %.1f dB\n', snr_db);
t = (0:numel(tx)-1)'/fs;
clean = h*tx(:).*exp(1i*2*pi*cfo_hz*t);
noise_sigma = sqrt(mean(abs(clean).^2)/(2*10^(snr_db/10)));
rx = clean + noise_sigma*(randn(size(clean))+1i*randn(size(clean)));

phyRx = LoRaPHY(rf,sf,bw,fs);
phyRx.fast_mode = true;
phyRx.cr = 4;
phyRx.crc = 1;
csi_opts = struct('guard_samples',2,'cfo_correction',true, ...
    'remove_static',true,'iter_opts',struct('pilot_fallback',true));
[csi,decoded,~,info] = phyRx.extract_csi(rx,5,csi_opts);

% Check payload bytes, then inspect confidence-selected CSI.
correct = ~isempty(decoded) && size(decoded,1)>=numel(payload) && ...
    isequal(uint8(decoded(1:numel(payload),1)),payload);
fprintf('Payload correctly decoded: %d\n',correct);
fprintf('Receiver packets detected: %d\n',numel(info));
if isempty(csi) || isempty(csi(1).corrected_csi)
    error('No usable confidence-symbol CSI from the received IQ.');
end
fprintf('Selected sensing chirps: %d\n',numel(csi(1).pilot_indices));
fprintf('Estimated packet CFO: %.3f Hz (injected %.1f Hz)\n', ...
    csi(1).cfo_est_hz,cfo_hz);

figure('Name','Integrated LoRa communication and CSI');
plot(csi(1).pilot_indices,abs(csi(1).corrected_csi),'o-','LineWidth',1.3);
grid on;
xlabel('Selected payload chirp index');
ylabel('Estimated |CSI|');
title(sprintf('Soft-decoded LoRa CSI (received SNR %.1f dB)', snr_db));
assert(correct,'The soft receiver did not recover the transmitted payload.');
