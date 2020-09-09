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


function APPEAR( EEG, config, outdir, suffix)

% Inputs 
%- config: Configuration input as structure 
%- subj_folder: subject folder that includes EEG files
%- subj_task_in: subject file to be cleaned. Should be located inside the scricpt  
% Load Input files


%% Define QA outputs 
% OZ: Define all constant here 
raw_first_ouput      = sprintf('%s/%s_p-0',outdir, suffix);

ECG_num              = config.ECG_ch_ind; % Channel # 32 is ECG recording 

QRS_filename         = sprintf('%s/%s_fMRIB_R.csv',outdir,...
    suffix); % output of QRS from 
QRS_fig             = sprintf('%s/%s_fMRIB_R.png',outdir,...
    suffix); % output of QRS from 
% From Chung KI code
R_filename           = sprintf('%s/%s_ICA_R.csv',outdir,...
    suffix);

R_fig                = sprintf('%s/%s_ICA_R.png',outdir,...
    suffix);% From Chung KI code
% Also output pulse ox QA
pulseocx_out_fig     = sprintf('%s/%s_pulse_oxc.png', outdir, suffix); % write unzipped file

pulseocx_output      = sprintf('%s/%s_pulse_oxc.csv',outdir,suffix); % write unzipped file

%OZ: ICA maps
corrEEG_ICA_maps     = sprintf('%s/%s_ICA_eeg_p-1.png',outdir, suffix);

%OZ: Final output file 
corrEEG_filename     = sprintf('%s/%s_eeg_p-1',outdir, suffix);

%OZ: Write averge HR calcuated from three methods 
HR_filename           = sprintf('%s/%s_HR_info.csv',outdir,...
    suffix);% From Chung KI code

% Also output ECG from EEG
ECG_EEG_filename      = sprintf('%s/%s_EEG_ECG.csv',outdir,...
    suffix);% From Chung KI code

polt_ecg_range         = 5:35 ;% Plot from 5 to 35 secs of ECG 


%% Gradient Artifact Correction, DownSampling, and Filtering
% DownSampling, Bandpass filtering, bandstop and sectioning the data
fprintf('Downsampling ..\n')
EEG = preproc(EEG);
% Downsmapling was set to 250 
% OZ: Save data after resmapling 
pop_writebva(EEG, strcat(outdir, '/', suffix, '_eeg_p-1'));


%% Try built-in BCG correction on EEGLab (FMRIB plugin)
% OZ: Let's clean unused varialbes and optimize the speed 

EEG1            = EEG;% to do QRS detection by fmrib
EEG2            = EEG;
% Apply some filtering 
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
% OZ: Calc. Heart Rate 
ECG                   = EEG1.data(ECG_num,:);
HR1                   = (Fs./(diff(QRS)))*60;
meanHR1               = mean(HR1);

% OZ: Write QRS file 
csvwrite(QRS_filename,QRS)


%OZ:  Also output ECG from EEG to csv
csvwrite(stract(outdir, '/', suffix, '_', 'fMRIB.csv'),ECG)

% Oz: Save figure of 30 secs
h=figure; plot(ECG)
hold on; plot( QRS ,...
    ECG(QRS),  'O');
xlim([polt_ecg_range(1)*EEG.APEAR.Fs,polt_ecg_range(end)*EEG.APEAR.Fs]);

saveas(h, stract(outdir, '/', suffix, '_', 'fMRIB.png'));
close(h);


%% Try CK QRS detection

% determine cardiac period
[pkloc,~, A, W, icaweights,icasphere] = cardperiod(EEG2,...
    EEG2.APEAR.slice_per_TR/EEG2.APEAR.TR,...
    outfile,outdir);

EEG2.event                            = AddMarker(EEG2.event,pkloc,1,0,...
    'R','Response');
EEG2.event                            = sortlatency(EEG2.event,1);

% OZ: Save Chung Ki sdetected Peaks 
csvwrite(R_filename,pkloc)

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
HR2           = (EEG2.APEAR.Fs./diff(Rpeaks))*60;
meanHR2       = mean(HR2);

h=figure; plot(ECG)
hold on; plot( Rpeaks,...
    ECG(Rpeaks),  'O');
xlim([polt_ecg_range(1)*Fs,polt_ecg_range(end)*Fs]);

saveas(h, stract(outdir, '/', suffix, '_', 'ChungKi.png'));
close(h);

%% Find the more accurate peak detection method
% OZ- cleaned code
if isfield(EEG.APEAR.PulseOX)
    pusleox_waveform = EEG.APEAR.PulseOX.pusleox_waveform;

    h=figure; plot(pusleox_waveform)
    hold on; plot( EEG.APEAR.PulseOX.Peaks,...
        pusleox_waveform(EEG.APEAR.PulseOX.Peaks),  'O');
    xlim([polt_ecg_range(1)*EEG.APEAR.PulseOX.Fs,...
        polt_ecg_range(end)*EEG.APEAR.PulseOX.Fs]);

    saveas(h, stract(outdir, '/', suffix, '_', 'PulseOx.png'));
    close(h);
    
    HR3          = (EEG.APEAR.Fs./diff(EEG.APEAR.PulseOX.Peaks))*60;% in min
    meanHR3      = mean(HR3);
    CK_MeanDif   = abs(meanHR3-meanHR2);
    % OZ: From Ahmad 
    fmrib_MeanDif= abs(meanHR3-meanHR1);
    
    %OZ: use pulse Ocx for correction 
    % Recommedned
    if usePulseocx
        EEG                    = fmrib_pas(EEG2,EEG.APEAR.PulseOX.Peak,'mean',21); %OZ: From Ahmad 
        selected               = 'PulseOX'; % For later analysis, mark that we selected Chung Ki ECG    
        
    else
        % Select either using Chung Ki's method or fMRIB 
        if  CK_MeanDif>(fmrib_MeanDif+0.1) %OZ:  From Ahmad
            % fMRIB
            EEG1     = pop_fmrib_pas(EEG1,'qrs','mean',21); %OZ:  From Ahmad
            %EEG3 = pop_fmrib_pas(EEG3,'qrs','obs',4);
            EEG      = EEG1;
            selected = 1; % For later analysis, mark that we selected fmrib ECG
            % 
        else
            % Chung Ki's method 
            EEG2     = fmrib_pas(EEG2,pkloc,'mean',21); %OZ: From Ahmad 
            %EEG3 = fmrib_pas(EEG3,pkloc,'obs',4); % From Ahmad 
            EEG      = EEG2;
            selected = 2; % For later analysis, mark that we selected Chung Ki ECG

        end
        
    end
    
    
else
    % Chung Ki's method 

    EEG2         = fmrib_pas(EEG2,pkloc,'mean',21); % Oz: from Ahamd
    %EEG3 = fmrib_pas(EEG3,pkloc,'obs',4);
    EEG          = EEG2;
    selected     = 2; % For later analysis, mark that we selected Chung Ki ECG

    
end

%OZ:  Write HR info to disck 
meanHR          = [meanHR1 meanHR2 meanHR3 selected];
writetable(stract(outdir, '/', suffix, '_', 'HR_info.png'), meanHR);

% eeg-ecg separation
[EEG.data,EEG.nbchan,~,~] = eegecgsep(EEG.data,EEG.nbchan,EEG.pnts,chlb);

clear EEG3 EEG1 EEG2
%% Bad Interval Detection
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

%% Run ICA
% OZ: I changed the code to run ICA once, since we already ran it in line 172 
n     = size(EEG.data,2);
    
EEG.pnts = n;
EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'interrupt','off');
EEG.chanlocs=loadbvef('BC-MR-32.bvef');
EEG.chanlocs(1)=[]; %removing GRND channel
EEG.chanlocs(1)=[]; %removing REF channel
EEG.chanlocs(32)=[]; %removing ECG channel
W = EEG.icaweights*EEG.icasphere;
A = inv(W);
EEG.icawinv     = pinv(EEG.icaweights*EEG.icasphere); % a priori same result as inv    

EEG.icachansind = 1:31;
EEG.chanlocs    = loadbvef('BC-MR-32.bvef');
EEG.chanlocs(1) = []; %removing GRND channel
EEG.chanlocs(1) = []; %removing REF channel
EEG.chanlocs(32)= []; %removing ECG channel

%OZ: plot IC maps
h              = figure;          
pop_topoplot(EEG, 0, [1:31] ,'',[6 6] ,0,'electrodes','on');
print(corrEEG_ICA_maps,'-dpng')

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
rate           = 1/Fs;
EEG.times      = [0:rate:(size(EEG.data,2)-1)*rate];

%% Classify ICs
%OZ: Obada cleaned the codes here 
% OZ: Ahmad set up the ICs to be chosen   
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

% Save final EEG
corrEEG_filename = stract(outdir, '/', suffix, '_', 'eeg_p-2')
% As edf
pop_writebva(finalEEG,corrEEG_filename);
% as MAT
save(strcat(corrEEG_filename, '.mat'),'finalEEG');
% As CSV
csvwrite(strcat(corrEEG_filename, '.csv'),finalEEG.data );

clear EEG_tmp; 
close all;
end







































