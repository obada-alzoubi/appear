% eeg-ecg signal separation
function [sgn,nbch,chlb,sgn1] = eegecgsep(sgn,nbch,smsz,chlb)
sgn1 = sgn(nbch,1:smsz);
sgn = sgn(1:nbch-1,1:smsz);
tmp = cell(nbch-1,3);
for ii=1:nbch-1
    tmp{ii,1} = chlb{ii,1};
end
chlb = tmp;
nbch = nbch-1;
end


