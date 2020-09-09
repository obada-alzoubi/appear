function [dst] = lne2pt(lnpt1, lnpt2, pt)
  A = lnpt2 - lnpt1;
  r = lnpt1 - pt;
  rxA = cross(r,A);
  magA = sqrt(sum(A.*A));
  dst = sqrt(sum(rxA.*rxA))/magA;
end
