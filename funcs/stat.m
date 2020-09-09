
% calculate few statistics
function [fave,fstd,fvar,fskw,fkurt] = stat(fn)
szfn = size(fn);
nr = szfn(1);
nc = szfn(2);
fave = zeros(1,nr);
fstd = zeros(1,nr);
fvar = zeros(1,nr);
fskw = zeros(1,nr);
fkurt = zeros(1,nr);
fave(:) = mean(fn(:,1:nc),2);
fstd(:) = std(fn(:,1:nc),0,2);
fvar(:) = var(fn(:,1:nc),0,2);
for ii=1:nr
    fskw(ii) = sum((fn(ii,1:nc)-fave(ii)).^3)/(nc*fvar(ii)*fstd(ii));
    fkurt(ii) = sum((fn(ii,1:nc)-fave(ii)).^4)/(nc*fvar(ii)*fvar(ii))-3;
end
end

