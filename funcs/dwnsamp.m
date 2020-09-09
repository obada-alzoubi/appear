
% downsampling
function dwndata = dwnsamp(data,orig,new)
dwndata = zeros(size(data,1)*orig/new,size(data,2));
for ii=1:size(data,2)
    dwndata(:,ii) = resample(data(:,ii),1,new/orig);
end
end
