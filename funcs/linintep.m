% bounded linear interpolation
function [yy0] = linintep(xx,yy,xx0)
yy0 = (yy(1)-yy(2))/(xx(1)-xx(2))*xx0 + (yy(2)*xx(1)-yy(1)*xx(2))/(xx(1)-xx(2));
if(yy0>max(yy))
    yy0 = max(yy);
elseif(yy0<min(yy))
    yy0 = min(yy);
end
end