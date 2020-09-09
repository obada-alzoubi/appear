%% This code is an attempt to automate all analysis from .vhdr up to IC classification
close all
clc
clear
%
% Information:
%      FMRIB plug-in: https://fsl.fmrib.ox.ac.uk/eeglab/fmribplugin/
%
%      Citation for icid.m: 
%      [Wong16] C.K. Wong, V. Zotev, M. Misaki, R. Philips, Q. Luo, and J. Bodurka,
%      (2016) Automatic EEG-assisted retrospective motion correction for fMRI (aE-REMCOR).
%      NeuroImage 129, pages 133-147.
% 
%      [Wong18] C.K. Wong, Q. Luo, V. Zotev, R. Phillips, K.W.C. Chan, and J. Bodurka,
%      (2018) Automatic cardiac cycle determination directly from EEG-fMRI data by 
%      multi-scale peak detection method.
%      submitted.
%    
%    
% Written by Kaylee Henry on 6/17/2019





% SET THESE BEFORE RUNNING THE CODE!!!!!
% 
% If you want the user to input:
% n=input('How many channels were used on the EEG?\n')
% ECG_chan=input('What channel was the ECG?\n')
% scntme=input('What was the scan time (in seconds)?\n')
% tr=input('What was the TR?\n')
% slpertr=input('How many slices per TR?\n')
% slmkpertr=input('How many markers are in each slice per TR?\n')
tic
n=32;
slpertr=42;
scntme=720;
tr=2;
slmkpertr=42;
ECG_chan=32;

indir='L:\jbodurka\Kaylee';
infile='AM395_ET32_25Sept180005';


% Step 1: Load in the EEG data as a .vhdr from Brain Analyzer
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    EEG = pop_loadbv(indir, 'AM395_ET32_25Sept180005.vhdr', [], [1:n]);
    
    % Finding the mrkA vector (using CK's function EEGLABB.m)
    [mrkA] = EEGLABB(indir,infile,scntme,tr,slmkpertr,slpertr);

         
% Step 2: MR gradient artifact removal & filtering
    EEG=preproc(EEG,tr,slpertr,mrkA);  
    
    
% Step 3: BCG Correction (FMRIB plugin)
    % Detect the QRS complexes
    EEG = pop_fmrib_qrsdetect(EEG,ECG_chan,'qrs','no');
    
    % Remove the BCG events
    EEG = pop_fmrib_pas(EEG,'qrs','obs',4);
    %eegplot(EEG.data)   
    
    
% Step 4: Raw data inspection (automatic continous cleaning)
    EEG = pop_rejcont(EEG, 'elecrange',[1:n-1] ,'freqlimit',[.5 70] ,'threshold',10,'epochlength',0.5,'contiguous',4,'addlength',0.25,'taper','hamming');
    %eegplot(EEG.data)
    %mixsig=double(EEG.data);
    
    
% Step 5: ICA
    % Run ICA
    EEG.data=EEG.data(1:n-1,:);
    EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'interrupt','on');
    %W=EEG.icaweights*EEG.icasphere; % seen in CK's EEGpreprocessing code
    %A=inv(W); % seen in CK's EEGpreprocessing code
    pop_eegplot(EEG,0,1,1) %plots the component
                           %pop_eegplot(EEG,1,1,1) plots the EEG
    
%    % Add the icachansind vector
    EEG.icachansind=1:n-1;
    
    % Load in the .locs file (make sure you change this so you have the correct number of channels loaded)
    EEG.chanlocs=loadbvef('BC-MR-32.bvef');
    EEG.chanlocs(1)=[]; %removing GRND channel
    EEG.chanlocs(1)=[]; %removing REF channel
    EEG.chanlocs(n)=[]; %removing ECG channel
    pop_topoplot(EEG, 0, [1:31] ,[6 6] ,0,'electrodes','on');
    %pop_prop( EEG, 0, 2, NaN, {'freqrange' [2 50] });
    
%% Step 6: ICLabel
    % Normal motion artifacts  
    EEG = iclabel(EEG, 'default');
    
    % The classification matrix is: EEG.etc.ic_classification.ICLabel.classifications
    [max_val, max_loc]=max(EEG.etc.ic_classification.ICLabel.classifications'); % The labels for the classification matrix are: EEG.etc.ic_classification.ICLabel.classes
    classmat=EEG.etc.ic_classification.ICLabel.classifications; %classification matrix
    % Analyzing the classifications
    lowclass=find(max_val<0.7); % find the max_values that are lower than 70%
    lowclass=lowclass(lowclass<33); % only considering channels #s lower than 33  
    for i=1:length(lowclass)
       if ((max_val(lowclass(1,i))<=.7) && (classmat(lowclass(1,i),1)>=0.15)) % if the classification is less than 70% and brain happens to be greater than 15%
           max_val(lowclass(1,i))=classmat(lowclass(i),1); %change the classification to be the brain
           max_loc(lowclass(1,i))=1;
       end
    end
    predictions=[max_loc; max_val]; % the values in a column labeled "1" are brain, anything else is an artifact
    toc
    
    % BCG IC identification (CK's function icid.m)
%     ic=double(EEG.data);
%     EEG.times=EEG.times/1000;
%     tic
%     cbicind=icid(ic,A,mixsig,EEG.srate,EEG.times(end));
%     toc %currently is taking 21 minutes to run icid()
  