% function max/min
function [pkl,lftmin,rgtmin,pkn,dpt,dptave,dptstd] = fninfo(fn)
[minv,minlc] = findpeaks(-fn);
nbmin = size(minv);
pkn = nbmin(2)-1;

pkl = zeros(1,pkn);
lftmin = zeros(1,pkn);
rgtmin = zeros(1,pkn);
lftmin(1:pkn) = minlc(1:pkn);
rgtmin(1:pkn) = minlc(2:pkn+1);
for ii = 1:pkn
    [~,pkl(ii)] = max(fn(lftmin(ii):rgtmin(ii)));
    pkl(ii) = pkl(ii) + lftmin(ii) - 1;
end
dpt = fn(pkl)-0.5*(fn(lftmin)+fn(rgtmin));

dptave = mean(dpt);
dptstd = std(dpt);
end


