% add new markers at the end of EEG.event
function NewEvent=AddMarker(event, MarkerTiming,duration,channel,type,code)

CurrEventSize=max(size(event,2),size(event,1));

for i=1:size(MarkerTiming,2)
    loc=CurrEventSize+i;
    event(loc).latency=MarkerTiming(i);
    event(loc).duration=duration;
    event(loc).channel=channel;
    event(loc).bvtime=[];
    event(loc).bvmknum=CurrEventSize+i;
    event(loc).type=type;
    event(loc).code=code;
    event(loc).urevent=CurrEventSize+i;
end
NewEvent=event;
end

