% bad segment definition
function badsegind = badseg(EEG,chlb)
[EEG.data,EEG.nbchan,~,~] = eegecgsep(EEG.data,EEG.nbchan,EEG.pnts,chlb);
ext = 0.3;
dscdt = 0.3;
sgn = EEG.data;
sgn = eegfilt(sgn,EEG.srate,1,0,0,3*fix(EEG.srate),0,'fir1',0);
nbch = EEG.nbchan;
smsz = EEG.pnts;
sgnmn = mean(EEG.data,2);
sgnstd = std(EEG.data,0,2);
durpt = floor(0.04*EEG.srate);
extpt = floor(ext*EEG.srate);
dscdpt = floor(dscdt*EEG.srate);
badsegind = zeros(size(sgn));
for ii=1:nbch
    ind = ( sgn(ii,:)>sgnmn(ii)+4*sgnstd(ii) | sgn(ii,:)<sgnmn(ii)-4*sgnstd(ii) );
    [sctn,sctr] = cont1(ind);
    if(sctn>0)
        for jj=1:sctn
            if(sctr(jj,2)-sctr(jj,1)>=durpt && max(abs(sgn(ii,sctr(jj,1):sctr(jj,2))))>200)
                sctr(jj,1) = max(1,sctr(jj,1) - extpt);
                sctr(jj,2) = min(smsz,sctr(jj,2) + extpt);
                ind(sctr(jj,1):sctr(jj,2)) = 1;
            else
                ind(sctr(jj,1):sctr(jj,2)) = 0;
            end
        end
        [sctn,sctr] = cont1(ind);
        if(sctn>0)
            for jj=1:sctn
                ltpt = sctr(jj,1);
                while( ltpt>1 && sgn(ii,ltpt)*sgn(ii,ltpt-1)>0 )
                    ind(ltpt-1) = 1;
                    ltpt = ltpt - 1;
                end
                rtpt = sctr(jj,2);
                while( rtpt<smsz && sgn(ii,rtpt)*sgn(ii,rtpt+1)>0 )
                    ind(rtpt+1) = 1;
                    rtpt = rtpt + 1;
                end
            end
            badsegind(ii,:) = ind;
        end
    end
end

for ii=1:nbch
    [sctn,sctr] = cont1(~badsegind(ii,:));
    for jj=1:sctn
        if(numel(sctr(jj,1):sctr(jj,2))<dscdpt)
            badsegind(ii,sctr(jj,1):sctr(jj,2)) = 1;
        end
    end
end
end

