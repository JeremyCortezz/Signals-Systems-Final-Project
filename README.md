# Signals-Systems-Final-Project
# Zot-Monitor-Filter

Wearable heart monitors are becoming increasingly important in biomedical engineering, allowing individuals to monitor cardiovascular health outside of clinical settings. Electrocardiograms (ECGs) provide critical insights into heart rate, rhythm, and overall cardiac activity. However, everyday physical activity introduces motion artifacts (MAs) caused by external body movement and internal physiological shifts, degrading signal quality and complicating analysis.

Traditional ECG noise-reduction techniques commonly rely on low-pass and high-pass filters to remove frequency-specific interference. While effective under relatively stationary conditions, these methods struggle during movement because motion artifacts can overlap with the frequency range of the underlying ECG signal. Consequently, filtering an entire recording uniformly can unintentionally distort important cardiac information.

The Zot Monitor Filter (ZMF) is a MATLAB-based signal-processing project designed to identify and reduce motion-related noise in wearable ECG recordings while preserving regions of the signal that remain relatively clean.

Project Objective

The goal of ZMF is to improve ECG signal quality during movement by detecting localized segments affected by motion artifacts and applying targeted noise-reduction techniques.

Rather than applying uniform filtering across an entire recording, ZMF aims to:

Identify regions of ECG data exhibiting elevated noise and motion artifacts.

Detect and characterize signal degradation using ECG and synchronized accelerometer data.

Apply noise reduction selectively to affected regions.

Preserve relatively clean ECG segments to minimize unnecessary signal distortion.

Improve the reliability of heart rate and cardiac rhythm analysis.

Methodology

The ZMF processing pipeline uses ECG data paired with synchronized accelerometer measurements to identify periods of motion-related interference.

The signal-processing workflow includes:

Signal preprocessing: Prepare ECG and accelerometer signals for analysis.

Artifact detection: Identify localized regions exhibiting characteristics associated with motion artifacts.

Noise characterization: Evaluate signal quality and distinguish relatively clean regions from noisy segments.

Selective filtering: Apply noise-reduction methods to detected artifact regions rather than uniformly filtering the entire recording.

Signal reconstruction: Combine processed regions with preserved portions of the original ECG.

Evaluation: Compare the original and filtered ECG signals to assess noise reduction and signal preservation.
