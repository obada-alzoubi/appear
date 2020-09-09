%% Script for loading in .edf ICA data from Brain Analyzer to validate the artifact labeling using EEGLAB's machine learning process

clc 
clear all
%close all

% Terminology:
    % EEG.icawinv is the inverse ICA weight matrix
    % EEG.icaweights is the ICA weight matrix
    % EEG.icachansind is a vector from 1:number of channels
    % EEG.icasphere is the sphere matrix created after pre-processing/whitening the data before ICA

% Loading in all of the files:
%files=dir('L:\jbodurka\EEG_Example\Outside_MR');

% Turn all of this into a for loop so we can run through multiple patients
    n=31; % the number of channels

% Load the .edf ICA data:
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    EEG = pop_biosig('L:\jbodurka\Kaylee\Raw_EEG\Processed\AP470_04\e\IC_Mat.edf', 'channels',[1:n]);
    indir = 'L:\jbodurka\Kaylee\Raw_EEG\Processed\AP470_04';
    rawdata = pop_loadbv(indir, 'AP470_AfterBCG.vhdr', [], [1:n]);
    EEG.times=EEG.times/1000;
    % AK762_05
    % AP470_04
    % AM359_08
    % AP552_04
    % AS098_06
    % AV134_02
    % AW686_05
    % AZ550_06
    % BC317_06
    % BC368_06
% Create the EEG.icachansind vector
    EEG.icachansind=1:n;

% Create an empty EEG.icasphere matrix
    EEG.icasphere=zeros(n,n);
 
% Load the ICA inverse matrix
    EEG.icawinv=load('L:\jbodurka\Kaylee\Raw_EEG\Processed\AP470_04\e\InvMixMat.txt');
    EEG.icawinv=EEG.icawinv(1:n,1:n);

% Load the ICA weights matrix
    EEG.icaweights=load('L:\jbodurka\Kaylee\Raw_EEG\Processed\AP470_04\e\MixMat.txt');
    EEG.icaweights=EEG.icaweights(1:n,1:n);
    
% Load in the .locs file
    %EEG=pop_chanedit(EEG, 'load',{'C:\\Users\\KHenry\\Desktop\\REU2\\eeglab2019_0\\BrainVision-10-20-Cap31.loc' 'filetype' 'autodetect'});
    EEG.chanlocs=loadbvef('BC-MR-32.bvef');
    EEG.chanlocs(1)=[]; %removing GRND channel?
    EEG.chanlocs(1)=[]; %removing REF channel?
    EEG.chanlocs(32)=[]; %removing ECG channel?
%
%
%
% Run ICLabel
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

    
% % % ICID
 EEG.icasphere=spher(EEG.data);
 W = EEG.icaweights*EEG.icasphere;
 A = inv(W);
 [cbicind,saccadey,blinky,topomap,spectrumy,tpblink,tpsac,smolregion,singchan,muscleloc] = icid(double(EEG.data),double(A),double(rawdata.data),EEG.srate,EEG.times(end));
 %uhh = makemap(double(W)*double(EEG.data),double(A));
