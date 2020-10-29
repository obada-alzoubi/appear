function [Peak_locations, pusleox_waveform] = pulseOx_DetectPeaks(pulseOx_file, pulseOx_Fs, EEG_Fs, minHearteRate)
    %% 
    % etect QRS indices for Pulse Ox. 
    % 
    % Read Pulse Ox. 
    pusleox_waveform     = csvread(pulseOx_file); % Read from pulse oxc 
    % Get sense of the data by estimating minimum peak height
    % Other ways for estimating minimum peaks height can be used.
    max1                 = max(pusleox_waveform(60*pulseOx_Fs:80*pulseOx_Fs));%
    max2                 = max(pusleox_waveform(200*pulseOx_Fs:220*pulseOx_Fs));%
    MinPeakAmp           = 0.25*min(max1,max2);  
    % Resamaple Pulse Ox to have same sampling rate as EEG   
    pusleox_waveform     = resample(pusleox_waveform,EEG_Fs, pulseOx_Fs ); 
    % Detect Peaks ( minimum heart rate to consider is minHearteRate)
    [~,Peak_locations]   = findpeaks(pusleox_waveform,'MinPeakDistance',...,
        round((minHearteRate*EEG_Fs)/60),'MinPeakHeight',MinPeakAmp);  
    
end

