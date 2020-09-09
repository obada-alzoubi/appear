function pow = power1(Fn,srt,tmwd,tmres,bndr)
fprintf('Calculating alpha power ...\n');
tm = tic;

wsz = 512;  ovlp = 256;  nfft = 1024;
PSD = zeros(size(Fn,1),nfft/2+1);
frq = zeros(1,nfft/2+1);
for ii=1:size(Fn,1)
    [PSD(ii,:),frq(:)] = pwelch(Fn(ii,:),wsz,ovlp,nfft,srt);
end
mpsd = mean(PSD,1);
[~,loc] = max(mpsd(frq>=bndr(1)&frq<=bndr(2)));
loc = loc + numel(frq(frq<bndr(1)));
pkfq = frq(loc);
pkrg = [pkfq-2,pkfq+2];

perlg = tmwd*srt;
stlg = tmres*srt;
sttt = floor((size(Fn,2)-perlg)/stlg);
pow = zeros(size(Fn,2)/(tmres*srt),size(Fn,1));
for ii = 1:sttt+1
    spt = zeros(perlg,size(Fn,1));
    frq = (1:1:perlg)/tmwd;
    for jj = 1:size(Fn,1)
        spt(:,jj) = abs(fft(transpose(Fn(jj,stlg*(ii-1)+1:stlg*(ii-1)+perlg).*transpose(hann(perlg))))/srt).^2;
    end
    pow(ii,:) = mean(spt(frq>=pkrg(1) & frq<=pkrg(2),:),1);
end
pow(sttt+2:end,:) = repmat(pow(sttt+1,:),size(pow(sttt+2:end,:),1),1);

fprintf('EEG alpha power calculation done -- ');
toc(tm);
fprintf('\n');
end

