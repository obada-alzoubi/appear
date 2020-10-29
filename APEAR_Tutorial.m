clear; clc;
% The toolbox was created at Laureate Institute for Brain Research 
% For technical issues, code and bugs please contact Obada Al Zoubi
% obada.y.alzoubi@gmail.com
% For sharing and general questions, please contact Dr. Jerzy Bodurka,
% Laureate Institute for Brain Research 
% This script is demo for running APPEAR Toolbox 

%% Add required paths and functions 
addpath('funcs');
addpath('eeglab2019_0');
[ALLEEG, ~, ~, ~] = eeglab;
close all;

sub_folder     =strcat(pwd, '/demo/');
% data are in brain vision analyzer format 
subj_eeg_file  = 'subj1_EEG';
% If you want to use pusle ox. to detect QRS peaks 
subj_ecg_file  = 'subj1_ECG.1D'; 

%% Output folder 

subj_out_folder = 'demo_out';
mkdir(subj_out_folder)

%% Store Configuration in EEG.APPEAR

TR              = 2;% seconds
slice_per_TR    = 39; % slices per volume 
scntme          = 480; % Scan length in secs

%% Read EEG, channel names and Slice Markers (R128) 

% Step 1
[EEG] = load_EEG(sub_folder,subj_eeg_file, scntme, TR, slice_per_TR);
% Set channel locations 
chanlocs                         = loadbvef('BC-MR-32.bvef');
chanlocs(1:2)                    = []; % Get rid of GND and REF
EEG.chanlocs                     = chanlocs;

EEG.APPEAR.Fs                    = 250;% Hz frequency of the EEG output 
EEG.APPEAR.filterRange           = [1 70]; %Hz - filter output EEG between 1 and 70 Hx
EEG.APPEAR.BCG_Crorrection       = 'Pulse_Ox'; % Recommended
EEG.APPEAR.PulseOx_Fs            = 40; % Hz
EEG.APPEAR.ECG_ch_ind            = 32; % ECG index (channel # 32)
EEG.APPEAR.PulseOX.minHearteRate = 25; % minimum heart rate you expect 
EEG.APPEAR.polt_ecg_range        = 5:35 ;% For QA, plot 30 sec of ECG waveform and the detected Peaks.

%% If we want to use Pulse Ox for BCG correction, add required fields  

% Step 2 
if ~isempty(subj_ecg_file) && strcmp(EEG.APPEAR.BCG_Crorrection, 'Pulse_Ox')
    % Read and resmaple Pulse Ox., then detect R peaks 
    [Peak_locations, pusleox_waveform]  = pulseOx_DetectPeaks(strcat(sub_folder, subj_ecg_file),...
      EEG.APPEAR.PulseOx_Fs  , EEG.APPEAR.Fs, EEG.APPEAR.PulseOX.minHearteRate) ;
   % Save QRS peaks, Fs and waveform  
    EEG.APPEAR.PulseOX.Peaks     = Peak_locations;
    EEG.APPEAR.PulseOX.Fs        = EEG.APPEAR.Fs;
    EEG.APPEAR.PulseOX.waveform  = pusleox_waveform ;
    
end

%% Run APPEAR 

EEG2 = APPEAR(EEG, subj_out_folder, 'test');

%% Output APPEAR to a file 

% Save final EEG
corrEEG_filename = strcat(subj_out_folder, '/', suffix, '_', 'eeg_p-2');
% As edf
pop_writebva(EEG2,corrEEG_filename);
% as MAT
save(strcat(corrEEG_filename, '.mat'),'EEG2');
% As CSV
csvwrite(strcat(corrEEG_filename, '.csv'),EEG2.data );