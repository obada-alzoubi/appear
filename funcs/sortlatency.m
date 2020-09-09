
% sort markers based on the duration timing 
function sorted = sortlatency(event,clmn2srt)
EventTable=struct2table(event);
SrtEventTable=sortrows(EventTable,clmn2srt);
NumEvent=1:1:size(SrtEventTable,1);
SrtEventTable.bvmknum=NumEvent';
SrtEventTable.urevent=NumEvent';
sorted=table2struct(SrtEventTable);
end

