clear;
nlink = 4;
psi = [.58541020,.13819660,.13819660,.13819660];
eta = [.13819660,.58541020,.13819660,.13819660];
chi = [.13819660,.13819660,.58541020,.13819660];
xe_sub = [1;1;0.5;1];
ye_sub = [1;1;0.5;0.5];
ze_sub = [1;0.5;0.5;0.5];
N1_sub = 0; N2_sub = 0; N3_sub = 0; N4_sub = 0;
fex = zeros(nlink,1);
xe = [1,1,1,0];
ye = [1,0,1,0];
ze = [1,1,0,0];
Amat = [ 1 1 1 1; xe(1) xe(2) xe(3) xe(4); ye(1) ye(2) ye(3) ye(4); ze(1) ze(2) ze(3) ze(4)];
volume_sub = tetraunitvolume('DUMMY')*parallelipipedvolume3d(xe_sub, ye_sub, ze_sub);

for quadpoint=1:4
    N1_sub = psi(quadpoint); N2_sub = eta(quadpoint);
    N3_sub = chi(quadpoint); N4_sub = 1-psi(quadpoint)-eta(quadpoint)-chi(quadpoint);

    xpos = N1_sub*xe_sub(1) + N2_sub*xe_sub(2) + N3_sub*xe_sub(3) + N4_sub*xe_sub(4);
    ypos = N1_sub*ye_sub(1) + N2_sub*ye_sub(2) + N3_sub*ye_sub(3) + N4_sub*ye_sub(4);
    zpos = N1_sub*ze_sub(1) + N2_sub*ze_sub(2) + N3_sub*ze_sub(3) + N4_sub*ze_sub(4);

    Npar = inv(Amat)*[1 xpos ypos zpos]';
%     uh = Npar(1)*u(1) + Npar(2)*u(2) + Npar(3)*u(3) + Npar(4)*u(4);
    %%%%%%% get exact solution %%%%%%%%%%%%%
%     [uex] = exactsolution(xpos,ypos,zpos);
    %%%%%%% get external force %%%%%%%%%%%%%
    fex = fex - 2*Npar*volume_sub*(1/4);
    %%%%%%% L2 Error Norm %%%%%%%%%%%%%%%%%
%     L2_error = L2_error + volume_sub*(uh-uex)*(uh-uex)*(1/4);
%     L2_exact = L2_exact + volume_sub*uex*uex*(1/4);
end