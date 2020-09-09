
function mbpfn = icsel(fcn,srate,cbind,acqfrq)
[nu,de] = butter(3,[2,2*acqfrq]/srate);
for ii=find(cbind==1)
    fcn(ii,:) = filtfilt(nu,de,fcn(ii,:));
end
sgn = zeros(1,size(fcn,1));
cnt = zeros(1,size(fcn,1));
lccnt = zeros(size(fcn,1),floor(size(fcn,2)/floor(4*srate)));
mn1 = mean(fcn,2);
std1 = std(fcn,[],2);
pssd = zeros(1,size(fcn,1));
ngsd = zeros(1,size(fcn,1));
for ii = find(cbind==1)
    fcn(ii,:) = fcn(ii,:)/std1(ii);
    mn1(ii) = mean(fcn(ii,:));
    std1(ii) = std(fcn(ii,:));
    pssd(ii) = mean( fcn(ii, fcn(ii,:)>=mn1(ii)+std1(ii) & fcn(ii,:)<=mn1(ii)+4*std1(ii) ) );
    ngsd(ii) = mean( fcn(ii, fcn(ii,:)<=mn1(ii)-std1(ii) & fcn(ii,:)>=mn1(ii)-4*std1(ii) ) );
    if(abs(pssd(ii))>=abs(ngsd(ii)))
        sgn(ii) = 1;
    elseif(abs(pssd(ii))<abs(ngsd(ii)))
        sgn(ii) = -1;
    end
    fcn(ii,:) = sgn(ii)*fcn(ii,:);
    mn1(ii) = sgn(ii)*mn1(ii);
    Fn = fcn(ii,:);  Fn(Fn>mn1(ii)+4*std1(ii))=mn1(ii)+4*std1(ii);
    for jj=1:floor(size(fcn,2)/floor(4*srate))
        seg = Fn((jj-1)*floor(4*srate)+1:jj*floor(4*srate));
        segpk = findpeaks(seg);
        segpk = sort(segpk,'descend');
        segval = sort(seg,'descend');
        lccnt(ii,jj) = (mean(segpk(1:round(0.1*numel(segpk))))-mean(segval(round(0.1*numel(segval)):end)))/std(segval(round(0.1*numel(segval)):end));
    end
    cnt(ii) = mean(lccnt(ii,:));
end
ind1 = find(cbind==1);
[~,ind] = max(cnt(ind1));
mbpfn = fcn(ind1(ind),:);
end
