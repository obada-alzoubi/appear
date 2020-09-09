%% Automating the EEG process up to ICLabel
clc
clear
close all
    n=31; % the number of channels

% Step 1: Load in the .vhdr data collected outside the scanner
% OR using .edf files
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    EEG = pop_biosig('L:\\jbodurka\\Kaylee\\RawData2.edf', 'channels',[1:n] );
    
% Step 2: Change the sampling rate from 5000 Hz to 250 Hz
    EEG = pop_resample(EEG, 250);
    
% Step 3: Bandpass filter [.5 80] Hz and 
    bpfrq = [.5,80];
    %fprintf('Band-pass frequency filtering the eeg data ...\n');
    bpfrq(2) = min(bpfrq(2),0.5*EEG.srate);
    EEG.data = eegfilt(EEG.data,EEG.srate,bpfrq(1),0,0,3*fix(EEG.srate/bpfrq(1)),0,'fir1',0);
    EEG.data = eegfilt(EEG.data,EEG.srate,0,bpfrq(2),0,3*fix(EEG.srate/bpfrq(1)),0,'fir1',0);

% Step 4: Raw data inspection using the automatic continous rejection feature in EEGLAB
    [EEG, selectedregions] = pop_rejcont(EEG, 'elecrange',[1:n] ,'freqlimit',[.5 80] ,'threshold',10,'epochlength',0.5,'contiguous',4,'addlength',0.25,'taper','hamming');
    
% Step 5: ICA
    EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'interrupt','on');
    %[wts,sph] = binica(EEG.data,'extended',1,'pca',31,'lrate',0.002,'stop',10^(-6),'maxsteps',1024,'verbose','off');
    %icasig = wts*sph*mixsig;
    %A = inv(wts*sph); 
% Step 6: Use ICtoICLabel.m to run ICLabel
    % Create the EEG.icachansind vector
    EEG.icachansind=1:n;
    
    % Load in the .locs file
     EEG.chanlocs=pop_chanedit(EEG, 'load',{'L:\\jbodurka\\Kaylee\\BrainVision-10-20-Cap31.loc' 'filetype' 'autodetect'});
     EEG.chanlocs=EEG.chanlocs.chanlocs;
%
%
    % Run ICLabel
     EEG = iclabel(EEG, 'default');

    % The classification matrix is: EEG.etc.ic_classification.ICLabel.classifications
     [max_val, max_loc]=max(EEG.etc.ic_classification.ICLabel.classifications'); % The labels for the classification matrix are: EEG.etc.ic_classification.ICLabel.classes
      predictions=[max_loc; max_val]; % the values in a column labeled "1" are brain, anything else is an artifact
      classmat=EEG.etc.ic_classification.ICLabel.classifications; %classification matrix
        
% Step 7: Analyzing the classifications
    lowclass=find(max_val<0.7); % find the max_values that are lower than 70%
    lowclass=lowclass(lowclass<33); % only considering channels #s lower than 33
    
    for i=1:length(lowclass)
       if ((max_val(lowclass(1,i))<=.7) && (classmat(lowclass(1,i),1)>=0.15)) % if the classification is less than 70% and brain happens to be greater than 15%
           max_val(lowclass(1,i))=classmat(lowclass(i),1); %change the classification to be the brain
       end
    end
    
% Step 8: Exporting the classifications & also the cleaned EEG ICA waveforms
    % NOTE: change the locations and names of the files you are saving!!
    % Save the weighted matrix as .txt file:
        %pop_expica(EEG, 'weights', 'L:\\jbodurka\\Kaylee\\weightICA');
        
    % Save the inverse weighted matrix as .txt file:
        %pop_expica(EEG, 'inv', 'L:\\jbodurka\\Kaylee\\invweightICA');
    
    % Export the data as .edf file:
        %pop_writeeeg(EEG, 'L:\\jbodurka\\Kaylee\\here', 'TYPE','EDF');
