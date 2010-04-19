function [L2_err_e, L2_ex_e, numintcheck, domvol] = L2error(node,j,B,x,y,z,disp,numnod,ifixu,volumep)

xe = x(node(1:4,j));
ye = y(node(1:4,j));
ze = z(node(1:4,j));
ue = disp(node(1:4,j));
Amat = [ 1 1 1 1; xe(1) xe(2) xe(3) xe(4); ye(1) ye(2) ye(3) ye(4); ze(1) ze(2) ze(3) ze(4)];
L2_err_e = 0; L2_ex_e = 0;

N1_sub = 0; N2_sub = 0; N3_sub = 0; N4_sub = 0;
N5_sub = 0; N6_sub = 0; N7_sub = 0; N8_sub = 0;
gaussfourpoint
% gauss = [-sqrt((3-2*sqrt(6/5))/7); sqrt((3-2*sqrt(6/5))/7); -sqrt((3+2*sqrt(6/5))/7); sqrt((3+2*sqrt(6/5))/7)];
% wt = [(18+sqrt(30))/36;(18+sqrt(30))/36; (18-sqrt(30))/36; (18-sqrt(30))/36];
% gauss = [-sqrt(3/5); 0; sqrt(3/5)];
% wt = [5/9;8/9;5/9];
voltet = 0;
numintcheck = 0;
for i=1:4
    for j=1:4
        for k=1:4

            psi = gauss(i); eta = gauss(j); chi = gauss(k);

            N1_sub = (1/8)*(1-psi)*(1-eta)*(1-chi);
            N2_sub = (1/8)*(1+psi)*(1-eta)*(1-chi);
            N3_sub = (1/8)*(1+psi)*(1+eta)*(1-chi);
            N4_sub = (1/8)*(1-psi)*(1+eta)*(1-chi);
            N5_sub = (1/8)*(1-psi)*(1-eta)*(1+chi);
            N6_sub = (1/8)*(1+psi)*(1-eta)*(1+chi);
            N7_sub = (1/8)*(1+psi)*(1+eta)*(1+chi);
            N8_sub = (1/8)*(1-psi)*(1+eta)*(1+chi);
            
            dNdpsi = 1/8*[-(1-eta)*(1-chi),(1-eta)*(1-chi),(1+eta)*(1-chi),-(1+eta)*(1-chi),-(1-eta)*(1+chi),(1-eta)*(1+chi),(1+eta)*(1+chi),-(1+eta)*(1+chi)];
            dNdeta = 1/8*[-(1-psi)*(1-chi),-(1+psi)*(1-chi),(1+psi)*(1-chi),(1-psi)*(1-chi),-(1-psi)*(1+chi),-(1+psi)*(1+chi),(1+psi)*(1+chi),(1-psi)*(1+chi)];
            dNdchi = 1/8*[-(1-eta)*(1-psi),-(1-eta)*(1+psi),-(1+eta)*(1+psi),-(1+eta)*(1-psi),(1-eta)*(1-psi),(1-eta)*(1+psi),(1+eta)*(1+psi),(1+eta)*(1-psi)];
            
            xesub = [xe(1);xe(1);xe(1);xe(1);xe(2);xe(3);xe(4);xe(2)];
            yesub = [ye(1);ye(1);ye(1);ye(1);ye(2);ye(3);ye(4);ye(2)];
            zesub = [ze(1);ze(1);ze(1);ze(1);ze(2);ze(3);ze(4);ze(2)];
            
            dxdpsi = dNdpsi*xesub; dxdeta = dNdeta*xesub; dxdchi = dNdchi*xesub;
            dydpsi = dNdpsi*yesub; dydeta = dNdeta*yesub; dydchi = dNdchi*yesub;
            dzdpsi = dNdpsi*zesub; dzdeta = dNdeta*zesub; dzdchi = dNdchi*zesub;
            
            J = [dxdpsi,dydpsi,dzdpsi
                dxdeta,dydeta,dzdeta
                dxdchi,dydchi,dzdchi];

            jcob_sub = abs(det(J));

            xpos = N1_sub*xe(1) + N2_sub*xe(1) + N3_sub*xe(1) + N4_sub*xe(1)...
                + N5_sub*xe(2) + N6_sub*xe(3) + N7_sub*xe(4) + N8_sub*xe(2);

            ypos = N1_sub*ye(1) + N2_sub*ye(1) + N3_sub*ye(1) + N4_sub*ye(1)...
                + N5_sub*ye(2) + N6_sub*ye(3) + N7_sub*ye(4) + N8_sub*ye(2);

            zpos = N1_sub*ze(1) + N2_sub*ze(1) + N3_sub*ze(1) + N4_sub*ze(1)...
                + N5_sub*ze(2) + N6_sub*ze(3) + N7_sub*ze(4) + N8_sub*ze(2);

            Npar = inv(Amat)*[1 xpos ypos zpos]';
            uh = Npar(1)*ue(1) + Npar(2)*ue(2) + Npar(3)*ue(3) + Npar(4)*ue(4);
            %%%%%%% Exact solution %%%%%%%%%%%%%
            [uex] = exactsolution(xpos,ypos,zpos);
            %%%%%%% L2 Error Norm %%%%%%%%%%%%%%%%%%
            L2_err_e = L2_err_e + wt(i)*wt(j)*wt(k)*jcob_sub*(uh-uex)*(uh-uex);
            L2_ex_e = L2_ex_e+ wt(i)*wt(j)*wt(k)*jcob_sub*uex*uex;
            voltet = voltet + wt(i)*wt(j)*wt(k)*jcob_sub;
            numintcheck = numintcheck + wt(i)*wt(j)*wt(k)*jcob_sub*uex;
        end
    end
end
domvol = voltet;