% Copyright (C) 2019 LaureateInstitute for Brain Research
% Change of 7/15/2020
% Code was edited by Obada Al Zoubi & Ahamd Mayeli, Laureate Institute for Brain Research
% Scripts were recorgnized to work with different platforms. 
% obada.y.alzoubi@gmail.com
% amayely@laureateinstitute.org
% Description: 
% This script runs EEG Preprocessing for EEG Recorded simultaneously with
% fMRI
% For preprocessing steps, please refer to Section 2.6:
% Mayeli, Ahmad, et al. "Real-time EEG artifact correction during fMRI using ICA." Journal of neuroscience methods 274 (2016): 27-37.
% It requires MATLAB and EEGLAB vs. eeglab2019_0
% Authors: Ahmad Mayeli, Kaylee Henry, Chung-Ki Wong
% Contributions: Dr. Jerzy Bodurka's Lab 
% http://www.laureateinstitute.org/jerzy-bodurka.html
% Citation for EEGLAB
%      Delorme A. & Makeig S. (2004), EEGLAB: an open source toolbox for analysis of single-trial EEG dynamics (pdf, 0.7 MB) Journal of Neuroscience Methods 134:9-21
%
% Citation for fmrib_qrsdetect.m and fmrib_pas.m
%   [Niazy06] R.K. Niazy, C.F. Beckmann, G.D. Iannetti, J.M. Brady, and
%   S.M. Smith (2005), Removal of FMRI environment artifacts from EEG data
%   using optimal basis sets. NeuroImage 28 (3), pages 720-737.
%   [Christov04] I.I. Christov (2004), Real time electrocardiogram QRS detection using 
%   combined adaptive threshold, BioMed. Eng. Online 3: 28.
%   [Kim04] K.H. Kim, et al., (2004), Improved ballistocardiac artifact removal 
%   from the electroencephalogram recored in fMRI, J NeouroSience Methods 135: 193-203.

%
% Please Cite:
%   Wong, Chung-Ki, et al. "Automatic cardiac cycle determination directly from EEG-fMRI data by multi-scale peak detection method." Journal of neuroscience methods 304 (2018): 168-184.
%   Mayeli, Ahmad, et al. "Real-time EEG artifact correction during fMRI using ICA." Journal of neuroscience methods 274 (2016): 27-37.
%   Wong, Chung-Ki, et al. "Automatic EEG-assisted retrospective motion correction for fMRI (aE-REMCOR)." Neuroimage 129 (2016): 133.


function [finalEEG] = APPEAR( EEG, outdir, suffix)

% Inputs 
%- EEG          : EEG object from EEGLAB with APPEAR structure 
%- outdir       : output directory
%- suffix       : suffix to add to the ouput files
%% General config. 

% plot only part of the ECG signal for QA. 
polt_ecg_range  = EEG.polt_ecg_range;% Plot from 5 to 35 secs of ECG 


%% Gradient Artifact Correction, DownSampling, and Filtering
% DownSampling, Bandpass filtering, bandstop and sectioning the data
fprintf('Downsampling ..\n')
EEG = preproc(EEG);
% Downsmapling was set to 250 
% OZ: Save data after resmapling 
pop_writebva(EEG, strcat(outdir, '/', suffix, '_eeg_p-1'));


%% Try built-in BCG correction on EEGLab (FMRIB plugin)
%Let's clean unused varialbes and optimize the speed 
% Make two copies of the EEG data to try different BCG corrections.
EEG1            = EEG;% to do QRS detection by fmrib
EEG2            = EEG;
% Apply some filtering for ECG 
ECG_num         = EEG.APPEAR.ECG_ch_ind; 
[nuECG,deECG]   = butter(3,[0.5,15]/(0.5*EEG.srate));
EEG1.data(ECG_num,:) = filtfilt(nuECG,deECG,EEG1.data(ECG_num,:));
EEG1            = pop_fmrib_qrsdetect(EEG1,ECG_num,'qrs','no'); % Detect the QRS complexes
evt             = EEG1.event;
evtsz           = size(evt);
mrkn            = 0;

%Calculate mean and std of heart rate using fmrib built in function for QRS
%Detection

for ii=1:max(evtsz(1),evtsz(2))
    
    if(strcmp(evt(ii).type,'qrs'))
        mrkn = mrkn + 1;        
        QRS(mrkn) = evt(ii).latency;
    end
    
end
% Calc. Heart Rate to compare different methods in detecting QRS 
ECG            = EEG1.data(ECG_num,:);
HR1            = (EEG1.srate./(diff(QRS)))*60;
meanHR1        = mean(HR1);

% Write QRS file fMRIB 
csvwrite(strcat(outdir, '/', suffix, '_', 'fMRIB_QRS.csv'),QRS)
% Also output ECG from EEG to csv
csvwrite(strcat(outdir, '/', suffix, '_', 'fMRIB_ECG.csv'),ECG)

% Save ECG figure for QA
h=figure; plot(ECG)
hold on; plot( QRS ,...
    ECG(QRS),  'O');
xlim([polt_ecg_range(1)*EEG.APPEAR.Fs,polt_ecg_range(end)*EEG.APPEAR.Fs]);

saveas(h, strcat(outdir, '/', suffix, '_', 'fMRIB.png'));
close(h);


%% Try Chung Ki's QRS detection
fprintf("Try to detect heart waveform peaks using Chung Ki's \n" )

% determine cardiac period
[pkloc,~, ~, ~, ~,~] = cardperiod(EEG2, EEG2.APPEAR.slice_per_TR/EEG2.APPEAR.TR, EEG.APPEAR.ECG_ch_ind);
% Addd R Peak markers from Chung Ki's methods
EEG2.event           = AddMarker(EEG2.event,pkloc,1,0, 'R','Response');
EEG2.event           = sortlatency(EEG2.event,1);

%Save Chung Ki sdetected Peaks for QA
csvwrite(strcat(outdir, '/', suffix, '_', 'ChungKi_Peaks.csv'),pkloc)

%Calculate mean and std of heart rate using CK function for QRS
%Detection
evt         = EEG2.event;
evtsz       = size(evt);
mrkn=0;
for ii=1:max(evtsz(1),evtsz(2))
    
    if(strcmp(evt(ii).type,'R'))
        mrkn = mrkn + 1;        
        Rpeaks(mrkn) = evt(ii).latency;
        
    end
end
% Calc. HR from Chung Ki method 
HR2           = (EEG2.APPEAR.Fs./diff(Rpeaks))*60;
meanHR2       = mean(HR2);
% Plot for QA
h=figure; plot(ECG)
hold on; plot( Rpeaks,...
    ECG(Rpeaks),  'O');
xlim([polt_ecg_range(1)*EEG2.APPEAR.Fs,polt_ecg_range(end)*EEG2.APPEAR.Fs]);

saveas(h, strcat(outdir, '/', suffix, '_', 'ChungKi.png'));
close(h);

%% Find  more accurate peak detection methods (use Pulse Ox.)

fprintf("Try to detect heart waveform peaks using pulse oximeter .. in case you have it \n" )
fprintf("Select best correction methods for BCG .. then apply it \n" )
% if you are sure that you want to use Pulse Ox. (recommended!)

switch EEG.APPEAR.BCG_Crorrection 
    
    case 'Pulse_Ox'
        pusleox_waveform = EEG.APPEAR.PulseOX.waveform;

        h=figure; plot(pusleox_waveform)
        hold on; plot( EEG.APPEAR.PulseOX.Peaks,...
            pusleox_waveform(EEG.APPEAR.PulseOX.Peaks),  'O');
        xlim([polt_ecg_range(1)*EEG.APPEAR.PulseOX.Fs,...
            polt_ecg_range(end)*EEG.APPEAR.PulseOX.Fs]);

        saveas(h, strcat(outdir, '/', suffix, '_', 'PulseOx.png'));
        close(h);

        HR         = (EEG.APPEAR.Fs./diff(EEG.APPEAR.PulseOX.Peaks))*60;% in min
        meanHR3     = mean(HR);
  
        EEG                    = fmrib_pas(EEG2,EEG.APPEAR.PulseOX.Peaks,'mean',21); %OZ: From Ahmad 
        selected               = 'PulseOX'; % For later analysis, mark that we selected Chung Ki ECG
    
    case 'fMRIB'         
           % fMRIB
            EEG     = pop_fmrib_pas(EEG1,'qrs','mean',21); %OZ:  From Ahmad
            %EEG3 = pop_fmrib_pas(EEG3,'qrs','obs',4);
            selected = 1; % For later analysis, mark that we selected fmrib ECG
            % 
    case 'MSPD'
            % Chung Ki's method 
            EEG     = fmrib_pas(EEG2,pkloc,'mean',21); %OZ: From Ahmad 
            %EEG3 = fmrib_pas(EEG3,pkloc,'obs',4); % From Ahmad 
            selected = 2; % For later analysis, mark that we selected Chung Ki ECG

    otherwise
        fprintf('Wrong method \n')
end 
%Write HR info to disck 
meanHR          = [meanHR1 meanHR2 meanHR3 selected];
writematrix(strcat(outdir, '/', suffix, '_', 'HR_info.png'), EEG.APPEAR.BCG_Crorrection);

% eeg-ecg separation
[EEG.data,EEG.nbchan,~,~] = eegecgsep(EEG.data,EEG.nbchan,EEG.pnts,EEG.APPEAR.chlb);

clear EEG3 EEG1 EEG2
%% Bad Interval Detection
fprintf("Deal with bad intervals and annotate them \n" )

% OZ: Cleaned the code 
EEG_tmp              = EEG; % save the EEG data to use later for matrix multiplication
SaveTime             = EEG.times;
EEG.xmax             = EEG.pnts/EEG.srate-1/EEG.srate;
EEG.chanlocs(ECG_num)= [];

% OZ: Ahmad set these Paramters
% OZ: pop_rejcont might not work 
badInter              =[];
[EEG_poten, badInter_poten]       = pop_rejcont(EEG, 'elecrange',[3:31] ,...
    'freqlimit',[0.5 7] ,'threshold',8,'epochlength',0.5,...
    'contiguous',4,'addlength',0.25,'taper','hamming');
% OZ: 
if ~isempty(EEG_poten.data)
    EEG     = EEG_poten;
    badInter= badInter_poten;
end

%% Run ICA for Examining ICs 
fprintf("Run ICA to classify  ICs \n" )

n     = size(EEG.data,2);
    
EEG.pnts = n;
EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'interrupt','off');

EEG.chanlocs(1)=[]; %removing GRND channel
EEG.chanlocs(1)=[]; %removing REF channel
EEG.chanlocs(EEG.APPEAR.ECG_ch_ind)=[]; %removing ECG channel
W = EEG.icaweights*EEG.icasphere;
A = inv(W);
EEG.icawinv     = pinv(EEG.icaweights*EEG.icasphere); % a priori same result as inv    

EEG.icachansind = setdiff(1:size(EEG.data,2), EEG.APPEAR.ECG_ch_ind);
EEG.chanlocs(1) = []; %removing GRND channel
EEG.chanlocs(1) = []; %removing REF channel
EEG.chanlocs(EEG.APPEAR.ECG_ch_ind)= []; %removing ECG channel

%OZ: plot IC maps
h              = figure;          
pop_topoplot(EEG, 0, setdiff(1:size(EEG.data,2), EEG.APPEAR.ECG_ch_ind) ,'',[6 6] ,0,'electrodes','on');
print(strcat(outdir, '/', suffix, '_', 'ICA_topo.png'),'-dpng')

close(h);

% x' = AS' where A is the mixing matrix and x' is the bad interval
% removed data, and S' is calculated below
S_prime = double(W)*double(EEG.data);

% Then we have x, A, S' and we want to solve x=AS for S
% So we do x*W = S
S               = double(W)*double(EEG_tmp.data);

% Now we need to get A' which is after we remove the artifactual IC
% OZ chnage from seconds to milliseconds 
EEG.times      = EEG.times/1000;
% OZ: Ahamd set these Parmaters
rate           = 1/EEG.APPEAR.Fs;
EEG.times      = [0:rate:(size(EEG.data,2)-1)*rate];

%% Classify ICs
%OZ: Obada cleaned the codes here 
% OZ: Ahmad set up the ICs to be chosen 
fprintf("Clean data based on the annotated ICs \n" )

[cbicind,~,~,~,~,tpblink,tpsac,...
    ~,singchan,muscleloc] = icid(W*double(EEG.data),...
    double(A),double(EEG.data),EEG.srate,EEG.times(end));

% Find the columns of A that have artifacts, and that gives us A'
tpblink   = logical(tpblink);
tpsac     = logical(tpsac);
muscleloc = logical(muscleloc);
singchan  = logical(singchan);

% blink or sac
bs        = or(tpblink,tpsac);
% muscle or sing chan
ms        = or(muscleloc,singchan');
% with bcg
allart    = bs+ms+cbicind;
allart    = sign(allart);
NeuralICs = find(allart==1);
% A to A'
A(:,NeuralICs) = 0;
A_prime   = A;
%% Reconstruct EEG with Inverse ICA
fprintf("Save cleaned data \n and exist \n" )

% Multiply A' by S, and that gives us the corrected EEG
EEG.times         = SaveTime;
finalEEG          = EEG;
finalEEG.times    = EEG_tmp.times;
finalEEG.data     = A_prime*S; %correct the times
finalEEG.pnts     = EEG_tmp.pnts;
finalEEG.event    = EEG_tmp.event;
finalEEG.chanlocs = EEG_tmp.chanlocs;

for i = 1:size(badInter,1)
    badlength      = badInter(i,2) - badInter(i,1);
    finalEEG.event = AddMarker(finalEEG.event,(badInter(i,1))',...
        badlength,0,'Userdefined','Bad Interval');
end

finalEEG.event     =sortlatency(finalEEG.event,1);

close all;


end







































