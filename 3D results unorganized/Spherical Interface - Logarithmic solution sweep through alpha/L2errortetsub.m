function [L2_error, L2_exact, domvol] = L2errortetsub(xe_sub,ye_sub,ze_sub,volume_sub,Amat,ue);

N1_sub = 0; N2_sub = 0; N3_sub = 0; N4_sub = 0;
N5_sub = 0; N6_sub = 0; N7_sub = 0; N8_sub = 0;
gaussfourpoint
% gauss = [-sqrt((3-2*sqrt(6/5))/7); sqrt((3-2*sqrt(6/5))/7); -sqrt((3+2*sqrt(6/5))/7); sqrt((3+2*sqrt(6/5))/7)];
% wt = [(18+sqrt(30))/36;(18+sqrt(30))/36; (18-sqrt(30))/36; (18-sqrt(30))/36];
% gauss = [-sqrt(3/5); 0; sqrt(3/5)];
% wt = [5/9;8/9;5/9];
L2_error = 0;
L2_exact = 0;
voltet = 0;
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
            
            xesub = [xe_sub(1);xe_sub(1);xe_sub(1);xe_sub(1);xe_sub(2);xe_sub(3);xe_sub(4);xe_sub(2)];
            yesub = [ye_sub(1);ye_sub(1);ye_sub(1);ye_sub(1);ye_sub(2);ye_sub(3);ye_sub(4);ye_sub(2)];
            zesub = [ze_sub(1);ze_sub(1);ze_sub(1);ze_sub(1);ze_sub(2);ze_sub(3);ze_sub(4);ze_sub(2)];
            
            dxdpsi = dNdpsi*xesub; dxdeta = dNdeta*xesub; dxdchi = dNdchi*xesub;
            dydpsi = dNdpsi*yesub; dydeta = dNdeta*yesub; dydchi = dNdchi*yesub;
            dzdpsi = dNdpsi*zesub; dzdeta = dNdeta*zesub; dzdchi = dNdchi*zesub;
            
            J = [dxdpsi,dydpsi,dzdpsi
                dxdeta,dydeta,dzdeta
                dxdchi,dydchi,dzdchi];

            jcob_sub = abs(det(J));

            xpos = N1_sub*xe_sub(1) + N2_sub*xe_sub(1) + N3_sub*xe_sub(1) + N4_sub*xe_sub(1)...
                + N5_sub*xe_sub(2) + N6_sub*xe_sub(3) + N7_sub*xe_sub(4) + N8_sub*xe_sub(2);

            ypos = N1_sub*ye_sub(1) + N2_sub*ye_sub(1) + N3_sub*ye_sub(1) + N4_sub*ye_sub(1)...
                + N5_sub*ye_sub(2) + N6_sub*ye_sub(3) + N7_sub*ye_sub(4) + N8_sub*ye_sub(2);

            zpos = N1_sub*ze_sub(1) + N2_sub*ze_sub(1) + N3_sub*ze_sub(1) + N4_sub*ze_sub(1)...
                + N5_sub*ze_sub(2) + N6_sub*ze_sub(3) + N7_sub*ze_sub(4) + N8_sub*ze_sub(2);

            Npar = inv(Amat)*[1 xpos ypos zpos]';
            uh = Npar(1)*ue(1) + Npar(2)*ue(2) + Npar(3)*ue(3) + Npar(4)*ue(4);
            %%%%%%% Exact solution %%%%%%%%%%%%%
            [uex] = exactsolution(xpos,ypos,zpos);
            %%%%%%% L2 Error Norm %%%%%%%%%%%%%%%%%%
            L2_error = L2_error + wt(i)*wt(j)*wt(k)*jcob_sub*(uh-uex)*(uh-uex);
            L2_exact = L2_exact + wt(i)*wt(j)*wt(k)*jcob_sub*uex*uex;
            voltet = voltet + wt(i)*wt(j)*wt(k)*jcob_sub;
        end
    end
end
domvol = voltet;