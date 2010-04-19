function [fex] = fextwedgesub(xe_sub,ye_sub,ze_sub,volume_sub,Amat);

N1_sub = 0; N2_sub = 0; N3_sub = 0; N4_sub = 0;
N5_sub = 0; N6_sub = 0; N7_sub = 0; N8_sub = 0;
% gauss = [-1/sqrt(3); 1/sqrt(3)];
gaussfourpoint
nlink = 4;
fex = zeros(nlink,1);
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
            
            xesub = [xe_sub(1);xe_sub(2);xe_sub(3);xe_sub(1);xe_sub(4);xe_sub(5);xe_sub(6);xe_sub(4)];
            yesub = [ye_sub(1);ye_sub(2);ye_sub(3);ye_sub(1);ye_sub(4);ye_sub(5);ye_sub(6);ye_sub(4)];
            zesub = [ze_sub(1);ze_sub(2);ze_sub(3);ze_sub(1);ze_sub(4);ze_sub(5);ze_sub(6);ze_sub(4)];
            
            dxdpsi = dNdpsi*xesub; dxdeta = dNdeta*xesub; dxdchi = dNdchi*xesub;
            dydpsi = dNdpsi*yesub; dydeta = dNdeta*yesub; dydchi = dNdchi*yesub;
            dzdpsi = dNdpsi*zesub; dzdeta = dNdeta*zesub; dzdchi = dNdchi*zesub;
            
            J = [dxdpsi,dydpsi,dzdpsi
                 dxdeta,dydeta,dzdeta
                 dxdchi,dydchi,dzdchi];
             
            jcob_sub = abs(det(J));

            xpos = N1_sub*xe_sub(1) + N2_sub*xe_sub(2) + N3_sub*xe_sub(3) + N4_sub*xe_sub(1)...
                + N5_sub*xe_sub(4) + N6_sub*xe_sub(5) + N7_sub*xe_sub(6) + N8_sub*xe_sub(4);

            ypos = N1_sub*ye_sub(1) + N2_sub*ye_sub(2) + N3_sub*ye_sub(3) + N4_sub*ye_sub(1)...
                + N5_sub*ye_sub(4) + N6_sub*ye_sub(5) + N7_sub*ye_sub(6) + N8_sub*ye_sub(4);

            zpos = N1_sub*ze_sub(1) + N2_sub*ze_sub(2) + N3_sub*ze_sub(3) + N4_sub*ze_sub(1)...
                + N5_sub*ze_sub(4) + N6_sub*ze_sub(5) + N7_sub*ze_sub(6) + N8_sub*ze_sub(4);
            
            rpos = xpos*xpos+ypos*ypos+zpos*zpos;

            Npar = inv(Amat)*[1 xpos ypos zpos]';
            %%%%%%% get External Force %%%%%%%%%%%
%             fex = fex - 6*wt(i)*wt(j)*wt(k)*zpos*Npar*jcob_sub;    %%%% Cubic solution
            %%%%%%% logarithmic solution %%%%%%%%%
            fex = fex - (1/rpos)*wt(i)*wt(j)*wt(k)*Npar*jcob_sub;
        end
    end
end