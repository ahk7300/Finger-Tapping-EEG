****OpenBCI EEG-Based Finger Tapping Interface for Stroke Rehabilitation***
This project explores an EEG-based Brain-Computer Interface (BCI) to decode motor intentions for aiding stroke rehabilitation. By analyzing alpha (8–12 Hz) and beta (13–30 Hz) waves, we distinguish between motor execution (ME) and motor imagery (MI) during finger tapping tasks. The project integrates signal processing and a gamified interface for enhanced neurorehabilitation.

Features
Real-time EEG signal processing with an 8-channel OpenBCI system.
Analysis of Event-Related Desynchronization (ERD) and Event-Related Synchronization (ERS).
TapTap Game interface for interactive stroke rehabilitation exercises.
Focus on motor execution and motor imagery for personalized recovery strategies.
Hardware and Setup
Equipment: OpenBCI Cyton board, EEG headset.
Electrode Placement:
Motor areas: C3, CZ, C4.
Parietal lobe: P3, PZ, P4.
Occipital lobe: O1, O2.
Data Collection:
Three conditions recorded for 3 minutes each:
Baseline (MB): No movement.
Motor Execution (ME): Actual finger tapping.
Motor Imagery (MI): Imagined finger tapping.
Signal Processing
Filters: Band-pass filter applied for alpha and beta bands.
Power Analysis: Spectral power computed using FFT for alpha and beta waves.
Key Metrics:
Alpha ERD: Decreased activity during ME/MI.
Beta ERS: Increased activity during motor execution.
TapTap Game
Interactive UI to classify ME and MI tasks.
Real-time feedback for stroke rehabilitation progress.
Threshold-based classification for motor recovery.
Results
Distinct ERD/ERS patterns observed for ME and MI tasks.
TapTap successfully differentiates between imagined and actual motor tasks.
Limitations
Connectivity issues with the OpenBCI system.
UI thresholds require refinement for better real-time accuracy.
Future Directions
Enhanced signal processing with adaptive filters.
Expanded UI for diverse motor impairments.
Broader applicability for neurorehabilitation beyond stroke recovery.
References
Refer to the detailed methodology and findings in the accompanying project report.

Authors
Yosemite Pinedo
Batool Elagha
Ali Haider Khan
Lorenzo Di Silverio
Danya Mohamed
