
# APPEAR Configuration Structure #
### APPEAR structure has the following fields that are needed for running the pipeline:
- APPEAR.Fs: A double, the desired downsampling for the output in Hz.
- APPEAR.FilterRange: A vector of two doubles, representing the lower and upper ranges to filter
the EEG data.
- APPEAR.BCG_Corrention: A string that takes either ‘Pulse_Ox’ or ‘fMRIb’ or ‘MSPD, the
desired method to detect QRS complex.
Pulse_Ox option uses pulse oximeter data to detect QRS complex locations. APPEAR also supports using
the fMRIb toolbox (Niazy et al., 2005) to detect QRS complex from the ECG electrode, if available.
Additionally, APPEAR supports detecting QRS complex from EEG data only if the ECG electrode was not
available, as in some cases. In this case, APPEAR uses a multi-scale peak detection method (MSPD) (Wong
et al., 2018). Please make sure you detect consistent and valid peaks when using any of the methods.
APPEAR will output the corresponding ECG signal with the detected peaks as quality assurance (QA).
- APPEAR.PulseOx_Fs: A double, the frequency of the pulse oximeter, if there is, and you wish to
use it for correction. In this case, APPEAR.BCG_Corrention should be set to ‘Pulse_Ox’
- APPEA.ECG_ch_ind: the index of the ECG channel (back electrode) in the EEG data if there is.
This is needed to inform APPEAR about the indices of EEG and ECG channels.
- APPEAR.PulseOx.minHearRate : A double, the minimum heart rate you expect when using
Pulse_Ox. It helps with detecting more accurate peaks.
If you are using ‘Pulse_Ox’, the following substructures will also be added.
- APPEAR.PulseOx.Peaks: A vector of integer, peaks location of BCG artifact
- APPEAR.PulseOx.Fs: A double, the sampling rate at which Pulse oximeter is used. We always
try to match the sampling rate of Pulse Oximeter with EEG frequency to ease the calculation.
This is done by either upsampling or downsampling Pulse Oximeter data.
- APPEAR.Pulse.Ox.waveform: A vector of doubles, the waveform of pulse oximeter after
resampling. APPEAR uses it for applying QA on the detected peaks. 
  - MRI scanner specific parameters:
  - TR: A double, time of repetition for the scanner.
  - slice_per_TR: An integer, the number of selected slices per TR.
  - scntme: An integer, scan time in secs.
  - Markers: An vector of integers, of Markers of each TR
- The previous information will be added to the APPEAR structure through load_EEG functions as
follows:
  - APPEAR.TR: A double, repletion time (e.g., 2)
  - APPEAR.slice_per_TR: An integer, slices per each TR (e.g., 39)
  - APPEAR.chlb: A cell of strings, channel labels
  - APPEAR.mrkA: A vector of integers, indices of each TR in the original sampling rate of EEG. It
should be equal to the number of TRs in the data. (e.g., TR=2 and scan time =200, then you
should have 100 values in the vector corresponding to the indices of the beginnings of each TR).
- APPEAR uses EEG object from EEGLAB to store APPEAR configuration. The structure will be under
EEG object with APPEAR name within EEG object:
```Matlab
>>EEG.APPEAR
``` 
<img src="images/p1.PNG" alt="hi" class="inline"/>
<img src="images/p2.PNG" alt="hi" class="inline"/>
