function [EEG,chlb,mrkA] = load_EEG(indir,infile,scntme,tr,slmkpertr)
% Inputs:
%    indir     - subject folder
%    infile    - subject file 
%    scntme    - length of the scan in secs
%    tr        - TR of the MR scanner
%    slmkpertr - Number of slices per TR
%
% Outputs:
%    EEG  - EEG object file
%    chlb - Channels label
%    mrkA - Markers of R128
%
% Author: Obada Al Zoubi
% Laureate Institute for Brain Research
% email:obada.y.alzoubi@gmail.com
% Website: http://www.obadaalzoubi.com
% Last revision: 7/19/2020

fprintf('Loading eeg data from %s/%s\n',indir,infile);
EEG = pop_loadbv(indir,strcat(infile,'.vhdr')); % loading EEG Data
% Extracting the EEG data features
srate = EEG.srate; %Sampling Rate
evt = EEG.event; %EEG Markers (Events)
evtsz = size(evt); % Number of Markers
chlb = chlin(); % Channel Names

% eeg-fmri data segmentation
smsz = scntme * srate;  % Number of EEG datapoints during scan
tic;
fprintf('Truncating eeg data ...\n');
mrkn = 0;
cs = floor(slmkpertr); % if number of slices per volume is not integer, find the diff
frs = slmkpertr - cs;
% find R128 markers and the last volume
for ii=evtsz(2):-1:1
    if(strcmp(evt(ii).type,'R128'))
        mrkn = mrkn + 1;
        cs = cs-1;
        if(cs==0)
            % tr2 is the begining of the last volume
            tr2 = evt(ii).latency - round(tr*srate*frs/slmkpertr);
        end
    end
end



tm2 = tr2 + tr*srate - 1;
tm1 = tm2 - smsz + 1;
mrkA = tm1:round(tr*EEG.srate):size(EEG.data,2);
mrkA = mrkA(1:round(scntme/tr));

EEG.bad = [];
EEG.badmot = [];
EEG.APEAR.TR              = tr;% seconds
EEG.APEAR.slice_per_TR    = slmkpertr; % slices per volume 
EEG.APEAR.scntme          = scntme; % Scan length in secs
EEG.APEAR.chlb            = chlb; % labels of channels 
EEG.APEAR.mrkA            = mrkA; % slice selection markers

toc;
end