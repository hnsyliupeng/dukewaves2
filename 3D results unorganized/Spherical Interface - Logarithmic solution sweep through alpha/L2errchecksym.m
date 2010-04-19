clear;
% readmesh
zdiv = 1
% alpha = 6*zdiv
% numtet
% r = 0.35;
[x,y,z,node,numtet,numnod] = manualmesh(zdiv);
% x = [1;1;1;0]; y = [1;0;1;0]; z = [2/3;2/3;1/3;1/3]; node=[1;2;3;4];
% numnod = 4; numtet = 1;
disp = zeros(numnod,1);
[ls] = calclevelset(x,y,z,numnod);
[ifixu, u] = essentialbcs(x,y,z,ls,numnod);
nlink=4;
% tol=1.0e-20;
% alpha = 30*zdiv;
keclass = zeros(nlink, nlink);
% initialize stiffness matrix and load vector
bigk = zeros(numnod,numnod);
fext = zeros(numnod,1);
ndof = 1;

for e=1:numtet
    %%%%% Calculate B Matrix %%%%%%%%%%
    [B] = sfderivatives(node,e,x,y,z);
    %%%%%%%%%%%% Test for intersection %%%%%%%%%%%%%%%%%%%%%%%%%
    [flag1,flag2] = intersectiontest(ls,node,e);
    if (flag1 == 1 && flag2 == 1)
%         fprintf('ls is zero')
        %%% Modify fixed nodes %%%%%%
        for nlink=1:4
            if(ls(node(nlink,e))<0)
                ifixu(node(nlink,e))=0;
                u(node(nlink,e)) = 0;
            end
        end
        %%%%% Calculate Intersection Points %%%%%%%%%
        [xint,yint,zint,coeff,nb_int,edge_id] = intersectionpointszerols(node,x,y,z,ls,e);
        %%%%% Calculate Normal %%%%%%%%%%%%%%%%%%%%%
        [normal_ls] = normal(coeff);
        %%%%% Calculate Volume on the positive side %%%%%%%%
        [volumep] = posvolumezerols(node,e,edge_id,nb_int,x,y,z,ls,xint,yint,zint);
        %%%%%%% Reorder intersection points %%%%%%%%%%%
        if(nb_int==4)
            [xint,yint,zint] = orderintpts(xint,yint,zint,normal_ls);
        end
        area_e = polyarea(xint,yint);
        if(volumep>0)
            alpha_e = (2*area_e)/(volumep);
        else
            alpha_e = 0;
        end
%         alpha(e) = alpha_e;
        %%%%% Calculate Nitsche and Penalty Terms %%%%%%%
        [kenit,fe_nit,kepen,fe_pen] = nitpenaltytermszerols(node,e,B,x,y,z,xint,yint,zint,nb_int,edge_id,normal_ls,volumep);
        
        %%%%% Calculate External Force Vector %%%%%%%%%%%%
        [fex] = force(node,e,B,x,y,z,disp,numnod,ls,ifixu,volumep);

    elseif(flag1==1 && flag2==0)
        for nlink=1:4
            if(ls(node(nlink,e))<0)
                ifixu(node(nlink,e))=0;
                u(node(nlink,e)) = 0;
            end
        end
        %%%%% Calculate Intersection Points %%%%%%%%%
        [xint,yint,zint,coeff,nb_int,edge_id] = intersectionpoints(node,x,y,z,ls,e);
        %%%%% Calculate Normal %%%%%%%%%%%%%%%%%%%%%
        [normal_ls] = normal(coeff);
        %%%%% Calculate Volume on the positive side %%%%%%%%
        [volumep] = posvolume(node,e,edge_id,nb_int,x,y,z,ls,xint,yint,zint);
        checkvolumep(e) = volumep;
        %%%%%%% Reorder intersection points %%%%%%%%%%%
        if(nb_int==4)
            [xint,yint,zint] = orderintpts(xint,yint,zint,normal_ls);
        end
        %%%%%%%%% Calculate Alpha %%%%%%%%%%%%%%%%%%%%%
        area_e = polyarea(xint,yint);
        if(volumep>0)
            alpha_e = (2*area_e)/(volumep);
        else
            alpha_e = 0;
        end
%         alpha(e) = alpha_e;
        %%%%% Calculate Nitsche and Penalty Terms %%%%%%%
        [kenit,fe_nit,kepen,fe_pen] = nitpenaltyterms(node,e,B,x,y,z,xint,yint,zint,nb_int,edge_id,normal_ls,volumep);
        %%%%% Calculate External Force Vector %%%%%%%%%%%
        [fex] = forcecutelem(node,e,B,x,y,z,xint,yint,zint,disp,numnod,ls,ifixu,nb_int,edge_id,volumep);

    else
        alpha_e = 0;
%         alpha(e) = alpha_e;
        %%%%%%%%%%%%% Volume of elements which are not intersected %%%%%%%%%%%%
        [volumep, kenit,fe_nit,kepen,fe_pen] = volcalc(node,e,x,y,z,ls);
        checkvolumep(e) = volumep;
        %%%%%%%%%%%%% External Force Vector %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [fex] = force(node,e,B,x,y,z,disp,numnod,ls,ifixu,volumep);
    end
    keclass = B'*B*volumep;
    ke = keclass + kenit + alpha_e*kepen;
    %assemble external force
    for i = 1:nlink
        fext(node(i,e)) = fext(node(i,e))+ fe_nit(i,1)...
            + alpha_e*fe_pen(i,1)+ fex(i,1);
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
[ifixu,u,nodevolumep] = nodeinfo(numtet,numnod,ls,node,x,y,z,ifixu,u);
bigk;
fext;
%%% Enforce Boundary conditions
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
K = sparse(bigk);
disp = K\fext;

%%%%%%%%%%% L2 error Symbolic %%%%%%%%%%%%%%%

syms u v t

xcoord = [1,1,1,0];
ycoord = [1,0,0,0];
zcoord = [0,1,0,0];

N = zeros(4,1);
A = [1,1,1,1;xcoord;ycoord;zcoord];
N = inv(A)*[1;1-t;1-t;t];
de = disp(node(1:4,1));
uh = N(1)*de(1)+N(2)*de(2)+N(3)*de(3)+N(4)*de(4);
L2errorsym_e = int(.5*(1-t)*(1-t)*(uh-(t*t+1))*(uh-(t*t+1)),t,0,1);