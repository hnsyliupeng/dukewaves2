function [xint,yint,zint] = orderthreeintpts(xint,yint,zint,normal_ls)

AB = [xint(2)-xint(1),yint(2)-yint(1),zint(2)-zint(1)];
BC = [xint(3)-xint(2),yint(3)-yint(2),zint(3)-zint(2)];
crossp1 = cross(AB,BC);

checknormal = dot(crossp1,normal_ls);

if (checknormal>0)
    xtemp = xint(2); ytemp = yint(2); ztemp = zint(2);
    xint(2) = xint(3); yint(2) = yint(3); zint(2) = zint(3);
    xint(3) = xtemp; yint(3) = ytemp; zint(3) = ztemp;
end