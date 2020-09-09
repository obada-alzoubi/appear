%% Script for loading in .edf ICA data from Brain Analyzer to validate the artifact labeling using EEGLAB's machine learning process

% Terminology:
    % EEG.icawinv is the inverse ICA weight matrix
    % EEG.icaweights is the ICA weight matrix
    % EEG.icachansind is a vector from 1:number of channels
    % EEG.icasphere is the sphere matrix created after pre-processing/whitening the data before ICA

% Loading in all of the files:
%files=dir('L:\jbodurka\EEG_Example\Outside_MR');

% Turn all of this into a for loop so we can run through multiple patients
    n=63; % the number of channels

% Load the .edf ICA data:
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    EEG = pop_biosig('L:\jbodurka\EEG_Example\Outside_MR\AI826\AI826_20160810\e\IC_Mat.edf', 'channels',[1:n]);

% Create the EEG.icachansind vector
    EEG.icachansind=1:n;

% Create an empty EEG.icasphere matrix
    EEG.icasphere=zeros(n,n);
 
% Load the ICA inverse matrix
    EEG.icawinv=load('L:\jbodurka\EEG_Example\Outside_MR\AI826\AI826_20160810\e\InvMixMat.txt');
    EEG.icawinv=EEG.icawinv(1:n,1:n);

% Load the ICA weights matrix
    EEG.icaweights=load('L:\jbodurka\EEG_Example\Outside_MR\AI826\AI826_20160810\e\MixMat.txt');
    EEG.icaweights=EEG.icaweights(1:n,1:n);
    
% Load in the .locs file
    %EEG=pop_chanedit(EEG, 'load',{'C:\\Users\\KHenry\\Desktop\\REU2\\eeglab2019_0\\BrainVision-10-20-Cap31.loc' 'filetype' 'autodetect'});
    EEG.chanlocs=loadbvef('BC-MR-64.bvef');
    EEG.chanlocs=EEG.chanlocs(1,3:65);
%
%
%
% Run ICLabel
    EEG = iclabel(EEG, 'default');

% The classification matrix is: EEG.etc.ic_classification.ICLabel.classifications
    [max_val, max_loc]=max(EEG.etc.ic_classification.ICLabel.classifications'); % The labels for the classification matrix are: EEG.etc.ic_classification.ICLabel.classes
    predictions=[max_loc; max_val]; % the values in a column labeled "1" are brain, anything else is an artifact
    classmat=EEG.etc.ic_classification.ICLabel.classifications; %classification matrix
    
%% Analyzing the classifications
    lowclass=find(max_val<0.7); % find the max_values that are lower than 70%
    lowclass=lowclass(lowclass<33); % only considering channels #s lower than 33
    
    for i=1:length(lowclass)
       if ((max_val(lowclass(1,i))<=.7) && (classmat(lowclass(1,i),1)>=0.15)) % if the classification is less than 70% and brain happens to be greater than 15%
           max_val(lowclass(1,i))=classmat(lowclass(i),1); %change the classification to be the brain
       end
    end