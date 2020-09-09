function [pkl,lftmin,rgtmin,pkn,dpt,dptave,dptstd] = fninfo(fn)
  [minv,minlc] = findpeaks(-fn);
  nbmin = size(minv);
  pkn = nbmin(2)-1;

  pkl = zeros(1,pkn);
  lftmin = zeros(1,pkn);
  rgtmin = zeros(1,pkn);
  lftmin(1:pkn) = minlc(1:pkn); %Vl
  rgtmin(1:pkn) = minlc(2:pkn+1); %Vr
  for ii = 1:pkn
     [~,pkl(ii)] = max(fn(lftmin(ii):rgtmin(ii))); %find max between two valleys
     pkl(ii) = pkl(ii) + lftmin(ii) - 1; %adjust the max location
  end
  dpt = fn(pkl)-0.5*(fn(lftmin)+fn(rgtmin)); % Vr<8 HzThe rise of the MO peaks as the average power diff between the peak and its neighboring right and left minima
%dpt is the rise of the MO for Vr<8Hz
  dptave = mean(dpt);
  dptstd = std(dpt);
end