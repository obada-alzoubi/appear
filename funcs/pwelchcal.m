function [psd,frq] = pwelchcal(fn,srate,wsz,ovl)
  nr = size(fn,1);
  nfft = 2*wsz;
  frq = zeros(1,nfft/2+1);
  psd = zeros(nr,nfft/2+1);
  for ii=1:nr
     [psd(ii,:),frq(:)] = pwelch(fn(ii,:),wsz,ovl,nfft,srate);
  end
end
