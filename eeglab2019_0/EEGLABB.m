function [mrkA] = EEGLABB(indir,infile,scntme,tr,slmkpertr,slpertr)

        % load eeg
          fprintf('Loading eeg data from %s/%s\n',indir,infile);
          EEG = pop_loadbv(indir,strcat(infile,'.vhdr'));
  
          srate = EEG.srate;
          evt = EEG.event;
          evtsz = size(evt); 
          

        % eeg-fmri data segmentation
          smsz = scntme * srate;
          tic;
          fprintf('Truncating eeg data ...\n');
          mrkn = 0;
          mrkA = zeros(evtsz(2));
          cs = floor(slmkpertr);
          frs = slmkpertr - cs;
          for ii=evtsz(2):-1:1
             if(strcmp(evt(ii).type,'R128'))
                mrkn = mrkn + 1;
                cs = cs-1;
                mrkA(ii) = evt(ii).latency;
                if(cs==0)
                   tr2 = evt(ii).latency - round(tr*srate*frs/slmkpertr);
                end
             end
          end
          mrkA = mrkA(mrkA~=0);
          tm2 = tr2 + tr*srate - 1;
          tm1 = tm2 - smsz + 1;
          mrkA = tm1:round(tr*EEG.srate):size(EEG.data,2);
          mrkA = mrkA(1:round(scntme/tr));
          toc;
end