% data segmentation
function [EEG,pkloc] = segmnt(EEG,pkloc,dscdt)
EEG.data = EEG.data(:,dscdt*EEG.srate+1:EEG.xmax);
EEG.bad = EEG.bad(:,dscdt*EEG.srate+1:EEG.xmax);
EEG.badmot = EEG.badmot(dscdt*EEG.srate+1:EEG.xmax);
pkloc = pkloc - round(dscdt*EEG.srate);
pkloc = pkloc(pkloc>0);
EEG.pnts = size(EEG.data,2);
EEG.xmax = size(EEG.data,2);
EEG.times = (1:1:EEG.xmax)/EEG.srate;
cnt=0;
for i=1:max(size(EEG.event,1),size(EEG.event,2))
    EEG.event(i).latency=EEG.event(i).latency-dscdt*EEG.srate;
    if EEG.event(i).latency<0
        cnt=cnt+1;
        DltMark(cnt)=i;
    end
end
if cnt>0
    EEG.event(DltMark)=[];
end
end