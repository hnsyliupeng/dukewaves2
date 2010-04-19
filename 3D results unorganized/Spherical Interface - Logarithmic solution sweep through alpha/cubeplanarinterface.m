%%%%%%%%% Solves Poisson's equation \del.(\del u) = 2 with a spherical
%%%%%%%%% inclusion of the form x*x + y*y + z*z = 0.45*0.45;
%%%%%%%%% Exact solution uex = z*z + 1
%%%%%%%%% Dirichlet boundaries z=1,y=1,x=1
%%%%%%%%% Neumann boundaries x=0, y=0, z=0 (Flux free)
clear;
area_int = pi*0.45*0.45*0.5;
volume_dom = 1-(1/6)*(pi*0.45*0.45*0.45);
C = area_int/volume_dom;
alpha_e = [C*(1e-02);C*(1e-01);2*C*(1e-01);3*C*(1e-01);4*C*(1e-01);5*C*(1e-01);6*C*(1e-01);7*C*(1e-01);8*C*(1e-01);9*C*(1e-01);C;2*C;3*C;4*C;5*C;6*C;7*C;
    8*C; 9*C; 10*C; 50*C; 70*C; C*(1e02); 2*C*(1e02); 3*C*(1e02); 5*C*(1e02); C*(1e03); C*(1e04); C*(1e5); C*(1e6);C*1e7; C*1e8; C*1e10; C*1e15; C*1e16];
for num_alpha = 1:size(alpha_e)
    readmesh
    zdiv = 4
    disp = zeros(numnod,1);
    [ls] = calclevelset(x,y,z,numnod);
    [ifixu, u] = essentialbcs(x,y,z,ls,numnod);
    nlink=4;
    % tol=1.0e-20;
    % alpha = 300*zdiv;
    keclass = zeros(nlink, nlink);
    % initialize stiffness matrix and load vector
    bigk = zeros(numnod);
    fext = zeros(numnod,1);
    ndof = 1;
    fprintf('Looping over elements \n')
    for e=1:numtet
        %%%%% Calculate B Matrix %%%%%%%%%%
        [B] = sfderivatives(node,e,x,y,z);
        %%%%%%%%%%%% Test for intersection %%%%%%%%%%%%%%%%%%%%%%%%%
        [flag1,flag2] = intersectiontest(ls,node,e);
        if (flag1 == 1 && flag2 == 1)
            %         fprintf('ls is zero')
            %         %%% Modify fixed nodes %%%%%%
            %         for nlink=1:4
            %             if(ls(node(nlink,e))<0)
            %                 ifixu(node(nlink,e))=0;
            %                 u(node(nlink,e)) = 0;
            %             end
            %         end
            %%%%% Calculate Intersection Points %%%%%%%%%
            [xint,yint,zint,coeff,nb_int,edge_id] = intersectionpointszerols(node,x,y,z,ls,e);
            %%%%% Calculate Normal %%%%%%%%%%%%%%%%%%%%%
            [normal_ls] = normal(coeff);
            %%%%% Calculate Volume on the positive side %%%%%%%%
            [volumep,volumeele] = posvolumezerols(node,e,edge_id,nb_int,x,y,z,ls,xint,yint,zint);
            %%%%%%% Reorder intersection points %%%%%%%%%%%
            if(nb_int==4)
                [xint,yint,zint] = orderintpts(xint,yint,zint,normal_ls);
            end
            area_e = polyarea(xint,yint);
%             alpha_e = 1*(area_int)/(volume_dom);
            %         alpha_e = 2*area_e/volumep;
            %         if(volumep>0)
            %             alpha_e = (2*area_e)/(volumep);
            %         else
            %             alpha_e = 0;
            %         end
            %         alpha(e) = alpha_e;
            %%%%% Calculate Nitsche and Penalty Terms %%%%%%%
            [kenit,fe_nit,kepen,fe_pen] = nitpenaltytermszerols(node,e,B,x,y,z,xint,yint,zint,nb_int,edge_id,normal_ls,volumep);

            %%%%% Calculate External Force Vector %%%%%%%%%%%%
            [fex] = force(node,e,B,x,y,z,disp,numnod,ls,ifixu,volumep);

        elseif(flag1==1 && flag2==0)
            %         for nlink=1:4
            %             if(ls(node(nlink,e))<0)
            %                 ifixu(node(nlink,e))=0;
            %                 u(node(nlink,e)) = 0;
            %             end
            %         end
            %%%%% Calculate Intersection Points %%%%%%%%%
            [xint,yint,zint,coeff,nb_int,edge_id] = intersectionpoints(node,x,y,z,ls,e);
            %%%%% Calculate Normal %%%%%%%%%%%%%%%%%%%%%
            [normal_ls] = normal(coeff);
            %%%%% Calculate Volume on the positive side %%%%%%%%
            [volumep,volumeele] = posvolume(node,e,edge_id,nb_int,x,y,z,ls,xint,yint,zint);
            checksmall(e) = volumep/volumeele;
            %%%%%%% Reorder intersection points %%%%%%%%%%%
            if(nb_int==4)
                [xint,yint,zint] = orderintpts(xint,yint,zint,normal_ls);
            end
            %%%%%%%%% Calculate Alpha %%%%%%%%%%%%%%%%%%%%%
            %         area_e = polyarea(xint,yint);
            %         if(volumep>0)
            %             alpha_e = (2*area_e)/(volumep);
            %         else
            %             alpha_e = 0;
            %         end
            %         alpha(e) = alpha_e;
            area_e = polyarea(xint,yint);
%             alpha_e = 1*(area_int)/(volume_dom);
            %         alpha_e = (2*area_e)/(volumep);
            %%%%% Calculate Nitsche and Penalty Terms %%%%%%%
            [kenit,fe_nit,kepen,fe_pen] = nitpenaltyterms(node,e,B,x,y,z,xint,yint,zint,nb_int,edge_id,normal_ls,volumep);
            %%%%% Calculate External Force Vector %%%%%%%%%%%
            e;
            [fex] = forcecutelem(node,e,B,x,y,z,xint,yint,zint,disp,numnod,ls,ifixu,nb_int,edge_id,volumep);

        else
            %         alpha_e = 0;
            %         alpha(e) = alpha_e;
            %%%%%%%%%%%%% Volume of elements which are not intersected %%%%%%%%%%%%
            [volumep, kenit,fe_nit,kepen,fe_pen] = volcalc(node,e,x,y,z,ls);
            checkvolumep(e) = volumep;
            %%%%%%%%%%%%% External Force Vector %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [fex] = force(node,e,B,x,y,z,disp,numnod,ls,ifixu,volumep);
        end
%         alpha(e) = alpha_e;
        keclass = B'*B*volumep;
        ke = keclass + kenit + alpha_e(num_alpha)*kepen;
        %     ke = keclass + alpha*kepen;
        %assemble external force
        for i = 1:nlink
            fext(node(i,e)) = fext(node(i,e))+ fe_nit(i,1)...
                + alpha_e(num_alpha)*fe_pen(i,1)+ fex(i,1);
            %         fext(node(i,e)) = fext(node(i,e))+ alpha*fe_pen(i,1)+ fex(i,1);

        end
        % assemble ke into bigk
        n1 = ndof-1;
        for i=1:nlink;
            for a=1:nlink;
                rbk = ndof*(node(i,e)-1) + 1;
                cbk = ndof*(node(a,e)-1) + 1;
                re = ndof*(i-1)+1;
                ce = ndof*(a-1)+1;
                bigk(rbk:rbk+n1, cbk:cbk+n1) = bigk(rbk:rbk+n1, cbk:cbk+n1) + ke(re:re+n1, ce:ce+n1);
            end
        end
    end
    fprintf('Eliminating nodes completely in the negative region \n')
    [ifixu,u,nodevolumep] = nodeinfo(numtet,numnod,ls,node,x,y,z,ifixu,u);
    %%% Enforce Boundary conditions
    fprintf('Enforcing boundary conditions \n')
    for n=1:numnod
        if (ifixu(n) == 1)
            for b=1:numnod
                fext(b) = fext(b) - bigk(b,n)*u(n);
            end
            bigk(n,1:numnod) = zeros(1,numnod);
            bigk(1:numnod,n) = zeros(1,numnod);
            bigk(n,n) = 1.0;
            fext(n) = u(n);
        end
    end
    [uexact] = exactsoln(x,y,z,numnod,ls,ifixu,u);
    % K = sparse(bigk);
    disp = bigk\fext;
    % [superr,err] = supnorm(disp, uexact, numnod,ifixu,x,y,z);
    % [superr,err] = supnormtest(disp,uexact,numnod,x,y,z,ls);
    % superr
    % % %%%%%%% If all subelements are tetrahedra %%%%%%%%%%%%%%%%%%%%%
    % L2err = 0;
    % L2ex = 0;
    % L2excheck = 0;
    % numintordercheck = 0;
    % for e=1:numtet
    %     xe = x(node(1:4,e));
    %     ye = y(node(1:4,e));
    %     ze = z(node(1:4,e));
    %     volumep = tetraunitvolume('DUMMY')*parallelipipedvolume3d(xe, ye, ze);
    % %     [volumep, kenit,fe_nit,kepen,fe_pen] = volcalc(node,e,x,y,z);
    %     [L2_err_e, L2_ex_e, numintcheck] = L2error(node,e,B,x,y,z,disp,numnod,ifixu,volumep);
    %     L2excheck = L2excheck + L2_ex_e;
    %     numintordercheck = numintordercheck + numintcheck;
    %     L2_err(e) = L2_err_e;
    %     L2err = L2err + L2_err_e;
    %     L2ex = L2ex + L2_ex_e;
    % end
    % L2norm = sqrt(L2err)/sqrt(L2ex)
    L2err = 0;
    L2ex = 0;
    L2excheck = 0;
    L2excheckcut = 0;
    numintordercheck = 0;
    numintelem = 0;
    domainvol = 0;
    domainvolcutelem = 0;
    domainvolnoncutelem = 0;
    volumecutelem = 0;
    fprintf('Calculating L2 norm \n')
    for e=numtet:-1:1
        originflag=0;
        [flag1,flag2] = intersectiontest(ls,node,e);
        if(flag1==1 && flag2==0)
            numintelem = numintelem + 1;
            [xint,yint,zint,coeff,nb_int,edge_id] = intersectionpoints(node,x,y,z,ls,e);
            [volumep] = posvolume(node,e,edge_id,nb_int,x,y,z,ls,xint,yint,zint);
            volumecutelem = volumecutelem + volumep;
            [normal_ls] = normal(coeff);
            if(nb_int==3)
                [xint,yint,zint] = orderthreeintpts(xint,yint,zint,normal_ls);
            elseif(nb_int==4)
                [xint,yint,zint] = orderintpts(xint,yint,zint,normal_ls);
            end
            e;
            [L2_err_e, L2_ex_e, domvol] = L2errorcutelem(node,e,B,x,y,z,xint,yint,zint,disp,numnod,ls,ifixu,nb_int,edge_id,volumep);
            L2excheckcut = L2excheckcut + L2_ex_e;
            domainvolcutelem = domainvolcutelem + domvol;
        elseif(flag1==0)
            [volumep, kenit,fe_nit,kepen,fe_pen] = volcalc(node,e,x,y,z,ls);
            if(volumep>0)
                [L2_err_e, L2_ex_e, numintcheck, domvol] = L2error(node,e,B,x,y,z,disp,numnod,ifixu,volumep);
                L2excheck = L2excheck + L2_ex_e;
                numintordercheck = numintordercheck + numintcheck;
                domainvolnoncutelem = domainvolnoncutelem + domvol;
            else
                L2_err_e = 0;
                L2_ex_e = 0;
            end
        end
        L2_err(e) = L2_err_e;
        L2err = L2err + L2_err_e;
        L2ex = L2ex + L2_ex_e;
    end
    L2norm(num_alpha) = sqrt(L2err)/sqrt(L2ex);
    % ExportField(x,y,z,node,disp)
end
plot(log(alpha_e),log(L2norm))