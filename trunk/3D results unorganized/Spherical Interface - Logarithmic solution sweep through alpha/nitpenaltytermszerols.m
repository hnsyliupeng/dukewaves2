function [kenit,fe_nit,kepen,fe_pen] = nitpenaltytermszerols(node,e,B,x,y,z,xint,yint,zint,nb_int,edge_id,normal_ls,volumep)

xe = x(node(1:4,e));
ye = y(node(1:4,e));
ze = z(node(1:4,e));
dNdx = B(1,1:4);
dNdy = B(2,1:4);
dNdz = B(3,1:4);

%2-pt gauss rule
gauss = [-3^(-0.5), 3^(-0.5)];
wtgp = [0.5, 0.5];
nlink = 4;
c1 = 1; c2 = 1;
kenit = zeros(4,4);
kepen = zeros(4,4);
fe_nit = zeros(nlink,1);
fe_pen = zeros(nlink,1);

Amat = [ 1 1 1 1; xe(1) xe(2) xe(3) xe(4); ye(1) ye(2) ye(3) ye(4); ze(1) ze(2) ze(3) ze(4)];

if(volumep<=0)
    kenit = zeros(nlink,nlink);
    kepen = zeros(nlink,nlink);
    fe_nit = zeros(nlink,1);
    fe_pen = zeros(nlink,1);
else

    if(nb_int==3)
        psi = [2/3 1/6 1/6]; eta = [1/6 2/3 1/6];

        for quadpoint=1:3

            N1_sub = 0; N2_sub = 0; N3_sub = 0; N4_sub = 0;

            dxdpsi = xint(1)-xint(3);
            dxdeta = xint(2)-xint(3);
            dydpsi = yint(1)-yint(3);
            dydeta = yint(2)-yint(3);
            dzdpsi = zint(1)-zint(3);
            dzdeta = zint(2)-zint(3);

            jcobxy = dxdpsi*dydeta - dydpsi*dxdeta;
            jcobyz = dydpsi*dzdeta - dydeta*dzdpsi;
            jcobxz = dzdpsi*dxdeta - dxdpsi*dzdeta;

            jcob_sub = sqrt(jcobxy*jcobxy + jcobxz*jcobxz + jcobyz*jcobyz);


            if (edge_id(1) == 21 & edge_id(2) == 31 & edge_id(3) == 41)

                N2_sub = psi(quadpoint);  N3_sub = eta(quadpoint); N4_sub = 1-psi(quadpoint)-eta(quadpoint);
                zpos = N2_sub*zint(1) + N3_sub*zint(2) + N4_sub*zint(3);
                ypos = N2_sub*yint(1) + N3_sub*yint(2) + N4_sub*yint(3);
                xpos = N2_sub*xint(1) + N3_sub*xint(2) + N4_sub*xint(3);

            elseif (edge_id(1) == 21 & edge_id(2) == 32 & edge_id(3) == 42)

                N1_sub = psi(quadpoint);  N3_sub = eta(quadpoint); N4_sub = 1-psi(quadpoint)-eta(quadpoint);
                zpos = N1_sub*zint(1) + N3_sub*zint(2) + N4_sub*zint(3);
                ypos = N1_sub*yint(1) + N3_sub*yint(2) + N4_sub*yint(3);
                xpos = N1_sub*xint(1) + N3_sub*xint(2) + N4_sub*xint(3);

            elseif (edge_id(1) == 31 & edge_id(2) == 32 & edge_id(3) == 43)

                N1_sub = psi(quadpoint);  N2_sub = eta(quadpoint); N4_sub = 1-psi(quadpoint)-eta(quadpoint);
                zpos = N1_sub*zint(1) + N2_sub*zint(2) + N4_sub*zint(3);
                ypos = N1_sub*yint(1) + N2_sub*yint(2) + N4_sub*yint(3);
                xpos = N1_sub*xint(1) + N2_sub*xint(2) + N4_sub*xint(3);

            elseif (edge_id(1) == 41 & edge_id(2) == 42 & edge_id(3) == 43)

                N1_sub = psi(quadpoint); N2_sub = eta(quadpoint); N3_sub = 1-psi(quadpoint)-eta(quadpoint);
                zpos = N1_sub*zint(1) + N2_sub*zint(2) + N3_sub*zint(3);
                ypos = N1_sub*yint(1) + N2_sub*yint(2) + N3_sub*yint(3);
                xpos = N1_sub*xint(1) + N2_sub*xint(2) + N3_sub*xint(3);
            end

            [uex] = exactsolution(xpos,ypos,zpos);
            Npar = inv(Amat)*[1 xpos ypos zpos]';

            dNdotN = dNdx*normal_ls(1) + dNdy*normal_ls(2) + dNdz*normal_ls(3);
            kenit = kenit - 0.5*(1/3)*jcob_sub*(Npar*dNdotN + dNdotN'*Npar');
            fe_nit = fe_nit - jcob_sub*(1/3)*(0.5)*(uex*dNdotN');

            %%%%%%%%%%%%%%%%%%%% Penalty Terms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            kepen = kepen + 0.5*(1/3)*jcob_sub*Npar*Npar';
            fe_pen = fe_pen + uex*jcob_sub*0.5*(1/3)*Npar;
        end

    elseif (nb_int == 4)
        
        for i=1:2
            for k=1:2
                psi = gauss(i); eta = gauss(k);

                N1_sub = 0.25*(1-psi)*(1-eta); N2_sub = 0.25*(1+psi)*(1-eta);
                N3_sub = 0.25*(1+psi)*(1+eta); N4_sub = 0.25*(1-psi)*(1+eta);
                dxdpsi = .25*(1-eta)*(xint(2)-xint(1)) + .25*(1+eta)*(xint(3)-xint(4));
                dxdeta = .25*(1-psi)*(xint(4)-xint(1)) + .25*(1+psi)*(xint(3)-xint(2));
                dydpsi = .25*(1-eta)*(yint(2)-yint(1)) + .25*(1+eta)*(yint(3)-yint(4));
                dydeta = .25*(1-psi)*(yint(4)-yint(1)) + .25*(1+psi)*(yint(3)-yint(2));
                dzdpsi = .25*(1-eta)*(zint(2)-zint(1)) + .25*(1+eta)*(zint(3)-zint(4));
                dzdeta = .25*(1-psi)*(zint(4)-zint(1)) + .25*(1+psi)*(zint(3)-zint(2));

                jcobxy = dxdpsi*dydeta - dydpsi*dxdeta;
                jcobxz = dzdpsi*dxdeta - dxdpsi*dzdeta;
                jcobyz = dydpsi*dzdeta - dzdpsi*dydeta;
                jcob_sub = sqrt( jcobxy*jcobxy + jcobxz*jcobxz + jcobyz*jcobyz);

                zpos = N1_sub*zint(1) + N2_sub*zint(2) + N3_sub*zint(3) + N4_sub*zint(4);
                ypos = N1_sub*yint(1) + N2_sub*yint(2) + N3_sub*yint(3) + N4_sub*yint(4);
                xpos = N1_sub*xint(1) + N2_sub*xint(2) + N3_sub*xint(3) + N4_sub*xint(4);

                [uex] = exactsolution(xpos,ypos,zpos);
                Npar = inv(Amat)*[1 xpos ypos zpos]';
                dNdotN = dNdx*normal_ls(1) + dNdy*normal_ls(2) + dNdz*normal_ls(3);
                kenit = kenit - jcob_sub*(Npar*dNdotN + dNdotN'*Npar');
                fe_nit = fe_nit - jcob_sub*(uex*dNdotN');

                %%%%%%%%%%%%%%%%%% Penalty terms %%%%%%%%%%%%%%%%%%%%%%%%%%
                kepen = kepen + jcob_sub*Npar*Npar';
                fe_pen = fe_pen + uex*jcob_sub*Npar;
            end
        end
    else
        kenit = zeros(nlink,nlink);
        fe_nit = zeros(nlink,1);
        kepen = zeros(nlink,nlink);
        fe_pen = zeros(nlink,1);
    end
end