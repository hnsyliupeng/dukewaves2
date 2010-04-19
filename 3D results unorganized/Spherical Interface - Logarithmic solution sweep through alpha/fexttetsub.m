function [fex] = fexttetsub(xe_sub,ye_sub,ze_sub,volume_sub,Amat);

nlink = 4;
psi = [.58541020,.13819660,.13819660,.13819660];
eta = [.13819660,.58541020,.13819660,.13819660];
chi = [.13819660,.13819660,.58541020,.13819660];
N1_sub = 0; N2_sub = 0; N3_sub = 0; N4_sub = 0;
fex = zeros(nlink,1);

for quadpoint=1:4
    N1_sub = psi(quadpoint); N2_sub = eta(quadpoint);
    N3_sub = chi(quadpoint); N4_sub = 1-psi(quadpoint)-eta(quadpoint)-chi(quadpoint);

    xpos = N1_sub*xe_sub(1) + N2_sub*xe_sub(2) + N3_sub*xe_sub(3) + N4_sub*xe_sub(4);
    ypos = N1_sub*ye_sub(1) + N2_sub*ye_sub(2) + N3_sub*ye_sub(3) + N4_sub*ye_sub(4);
    zpos = N1_sub*ze_sub(1) + N2_sub*ze_sub(2) + N3_sub*ze_sub(3) + N4_sub*ze_sub(4);

    Npar = inv(Amat)*[1 xpos ypos zpos]';
    
    rpos = xpos*xpos+ypos*ypos+zpos*zpos;
%     fex = fex - 6*zpos*Npar*volume_sub*(1/4);   %%%% Cubic solution
    fex = fex - (1/rpos)*Npar*volume_sub*(1/4); %%%% Logarithmic solution
end