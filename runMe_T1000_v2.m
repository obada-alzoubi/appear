clear; clc;
% The toolbox was created at Laureate Institute for Brain Research 
% For technicalissues, code and bugs please contact Obada Al Zoubi
% obada.y.alzoubi@gmail.com
% For sharing and general questions, please contact Dr. Jerzy Bodurka,
% Laureate Institute for Brain Research 
% This script is demo for running APEAR method 

%% Add required paths and functions 
addpath('funcs');
addpath('eeglab2019_0');
[ALLEEG, ~, ~, ~] = eeglab;
close all;


%% Input data 

sub_folder     =strcat(pwd, '/demo/');
% data are in brain vision analyzer format 
subj_eeg_file  = 'subj1_EEG';
% If you want to use pusle ox. to detect QRS peaks 
subj_ecg_file  = 'subj1_ECG.1D'; 

%% Output folder 

subj_out_folder = 'demo_out';
mkdir(subj_out_folder)

%% Store Configuration in EEG.APEAR
TR              = 2;% seconds
slice_per_TR    = 39; % slices per volume 
scntme          = 480; % Scan length in secs

%% Read EEG, channel names and Slice Markers (R128) 
% Step 1
[EEG] = load_EEG(sub_folder,subj_eeg_file, scntme, TR, slice_per_TR);

EEG.APEAR.Fs                   = 250;% Hz frequency of the EEG oupt 
EEG.APEAR.filterRange           = [1 70]; %Hz - filter output EEG between 1 and 70 Hx
EEG.APEAR.BCG_Crorrection       = 'Pulse_Ox'; % Recommended
EEG.APEAR.PulseOx_Fs            = 40; % Hz
EEG.APEAR.ECG_ch_ind            = 32; % ECG index (channel # 32)
EEG.APEAR.PulseOX.minHearteRate = 25; % minimum heart rate you expect 

%% If we want to use Pulse Ox for BCG correction, add required fields  
% Step 2 
if ~isempty(subj_ecg_file) && strcmp(EEG.APEAR.BCG_Crorrection, 'Pulse_Ox')
    % Read and resmaple Pulse Ox., then detect R peaks 
    [Peak_locations, pusleox_waveform]  = pulseOx_DetectPeaks(strcat(sub_folder, subj_ecg_file),...
      EEG.APEAR.PulseOx_Fs  , EEG.APEAR.Fs, EEG.APEAR.PulseOX.minHearteRate) ;
   % Save QRS peaks, Fs and waveform  
    EEG.APEAR.PulseOX.Peaks     = Peak_locations;
    EEG.APEAR.PulseOX.Fs        = EEG.APEAR.Fs;
    EEG.APEAR.PulseOX.waveform  = pusleox_waveform ;
    
end


%% Run APEAR 
status = APEAR(EEG, subj_out_folder, 'test');

