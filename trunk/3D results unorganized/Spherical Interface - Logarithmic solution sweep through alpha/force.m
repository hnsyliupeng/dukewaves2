function [fex] = force(node,e,B,x,y,z,disp,numnod,ls,ifixu,volumep)

xe = x(node(1:4,e));
ye = y(node(1:4,e));
ze = z(node(1:4,e));
u = disp(node(1:4,e));
nlink = 4;

psi = [.58541020,.13819660,.13819660,.13819660];
eta = [.13819660,.58541020,.13819660,.13819660];
chi = [.13819660,.13819660,.58541020,.13819660];
volume_sub = volumep;
fex = zeros(nlink,1);

N1_sub = 0; N2_sub = 0; N3_sub = 0; N4_sub = 0;
Amat = [ 1 1 1 1; xe(1) xe(2) xe(3) xe(4); ye(1) ye(2) ye(3) ye(4); ze(1) ze(2) ze(3) ze(4)];
for quadpoint=1:4

    N1_sub = psi(quadpoint); N2_sub = eta(quadpoint);
    N3_sub = chi(quadpoint); N4_sub = 1-psi(quadpoint)-eta(quadpoint)-chi(quadpoint);

    xpos = N1_sub*xe(1) + N2_sub*xe(2) + N3_sub*xe(3) + N4_sub*xe(4);
    ypos = N1_sub*ye(1) + N2_sub*ye(2) + N3_sub*ye(3) + N4_sub*ye(4);
    zpos = N1_sub*ze(1) + N2_sub*ze(2) + N3_sub*ze(3) + N4_sub*ze(4);
    
    rpos = xpos*xpos+ypos*ypos+zpos*zpos;
    Npar = inv(Amat)*[1 xpos ypos zpos]';
    %%%%%%%%% External Force Vector %%%%%%%%%%%%%%%%%%%%%%%%%%%
%     fex = fex - 6*zpos*Npar*volume_sub*(1/4);     %%%% Cubic solution
    fex = fex - (1/rpos)*Npar*volume_sub*(1/4);   %%%% Logarithmic solution
    
%     uh = Npar(1)*u(1) + Npar(2)*u(2) + Npar(3)*u(3) + Npar(4)*u(4);
    %%%%%%% get exact solution %%%%%%%%%%%%%
%     [uex] = exactsolution(xpos,ypos,zpos);
    %%%%%%% L2 Error Norm %%%%%%%%%%%%%%%%%
%     L2_error = L2_error + volume_sub*(uh-uex)*(uh-uex)*(1/4);
%     L2_exact = L2_exact + volume_sub*uex*uex*(1/4);
end