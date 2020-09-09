function [myPsd] = powerSpectral(X)
h1=spectrum.welch;
set(h1,'Windowname','Hann');
Fs=250;
set(h1,'OverlapPercent',66.7);
set(h1,'SegmentLength',512);
myPsd=psd(h1,X-mean(X),'Fs',Fs);


end

