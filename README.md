# SenLoRa-inspired MATLAB — essential files

A compact MATLAB implementation of integrated LoRa **communication and sensing**, inspired by the SenLoRa paper. The package is **files only**: the five core `.m` files and practical examples. 

**This is an independent research implementation, not official SenLoRa code or a drop-in LoRaWAN receiver.** It is built on the original LoRaPHY implementation by jkadbear. The original hard-only decoding entry points have been removed; the remaining shared synchronization/PHY helpers are still needed by the soft receiver.

## Files

| File | Purpose |
| --- | --- |
| `LoRaPHY.m` | LoRa modulation, synchronization, iterative soft decoding, confidence-symbol selection, complex CSI extraction, wrappers for activity/entropy/incentive/closed-loop processing |
| `senlora_activity.m` | Walking score and respiration-rate fitting from packet-timestamped CSI amplitude |
| `senlora_sample_entropy.m` | Sample entropy from the SenLoRa appendix, with matching-pair diagnostics |
| `senlora_incentive.m` | Offline incentive-decision analysis |
| `senlora_closed_loop.m` | Generate and process actual simulated IQ from periodic and requested extra packets |

### Examples

- `examples/communication_and_csi.m` — transmit a payload, recover it using the soft decoder, estimate CFO and plot confidence-chirp CSI.
- `examples/breathing_18bpm.m` — **validated** 18-BPM demonstration using real simulated LoRa IQ, matched ground truth and a breathing-frequency fit. User previously reported 46/46 payloads, 46/46 usable CSI packets, 18.000 BPM recovered at 20 dB received SNR.
- `examples/walking_vs_quiet.m` — compare two complete IQ receiver runs: quiet channel versus a moving reflector.
- `examples/adaptive_packet_requests.m` — closed-loop simulation that actually creates and decodes requested extra IQ packets when the controller triggers. Requests can legitimately be zero.

## Requirements

- MATLAB **R2019b or newer**.


## Running

Open MATLAB and set **Current Folder** to the directory containing `LoRaPHY.m` and this README.

```matlab
clear classes
which LoRaPHY -all

run('examples/communication_and_csi.m')
run('examples/breathing_18bpm.m')
run('examples/walking_vs_quiet.m')
run('examples/adaptive_packet_requests.m')
```

The examples are independent; run any one of them. The 18-BPM example takes longer than the single-packet example because it processes 46 packets.

## What is implemented

- **Communication:** original LoRa waveform/encoding helpers plus soft demodulation and Hamming decoding, stabilized log-sum-exp bit LLRs, and iterative confidence-based refinement.
- **Sensing:** confidence-symbol reference-signal reconstruction and complex CSI extraction; packet-mean CSI amplitude, motion classification and respiration-frequency fitting.
- **Adaptive operation:** sample entropy and simulated packet requests, with causal processing of newly generated IQ packets.

## Important research limitations

- Some `amendLLR`/channel-refinement details are **independent approximations**: the exact authors' Algorithm 2 implementation was not available. An alternative Python posterior decoder from previous experiments is **not silently substituted** for the tested MATLAB receiver.
- Phase 3 does not fully reproduce joint fine CFO/STO estimation, hardware SFO correction or PLL phase-jump interpolation from the paper. It estimates packet-relative CFO for simulation.
- `pilot_fallback=true` in the closed-loop simulator is an explicitly labeled **simulation heuristic** and can add low-confidence pilots.
- Sample entropy uses `m=2`, `r=0.2 × std` by default. **The paper did not publish these numerical settings**. Its experimental `gamma=0.01496` is a reference, not a generally calibrated threshold.
- Motion thresholds, channel paths, packet cadence and noise are simulation assumptions. The incentive code sends **simulated IQ**, not actual LoRaWAN signaling. Do not claim reproduction of the authors' hardware evaluation from these files.

## Attribution and redistribution

The original `LoRaPHY.m` was published by **jkadbear**, with an MIT License notice, at <https://github.com/jkadbear/LoRaPHY>, corresponding to *From Demodulation to Decoding: Toward Complete LoRa PHY Understanding and Implementation* (Xu et al.). See `THIRD_PARTY_LICENSE.txt` for the original code's MIT notice.

This package is **inspired by** *SenLoRa: Integrated Sensing and Communication with Ambient LoRa* (Wang et al., Proc. ACM IMWUT, 2025). It is not affiliated with the SenLoRa authors. The license terms for any new contribution must be decided by its copyright owner before publication; the original LoRaPHY MIT license must be preserved.
