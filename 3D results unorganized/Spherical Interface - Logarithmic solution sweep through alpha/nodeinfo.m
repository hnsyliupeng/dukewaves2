function [ifixu,u,nodevolumep] = nodeinfo(numtet,numnod,ls,node,x,y,z,ifixu,u)

nodevolumep = zeros(numnod,1);
tol = 1.0e-06;
nlink=4;
for e=1:numtet
    [flag1,flag2] = intersectiontest(ls,node,e);
    if(flag1==1 && flag2==1)
        [xint,yint,zint,coeff,nb_int,edge_id] = intersectionpointszerols(node,x,y,z,ls,e);
        [volumep]=posvolumezerols(node,e,edge_id,nb_int,x,y,z,ls,xint,yint,zint);
        for i=1:nlink
            nodevolumep(node(i,e)) = nodevolumep(node(i,e)) + volumep;
        end
    elseif(flag1==1 && flag2==0)
        [xint,yint,zint,coeff,nb_int,edge_id] = intersectionpoints(node,x,y,z,ls,e);
        [volumep]= posvolume(node,e,edge_id,nb_int,x,y,z,ls,xint,yint,zint);
        for i=1:nlink
            nodevolumep(node(i,e)) = nodevolumep(node(i,e)) + volumep;
        end
    elseif(flag1==0)
        [volumep, kenit,fe_nit,kepen,fe_pen] = volcalc(node,e,x,y,z,ls);
        for i=1:nlink
            nodevolumep(node(i,e)) = nodevolumep(node(i,e)) + volumep;
        end
    end
end

for nod=1:numnod
    if(nodevolumep(nod)>0)
        ifixu(nod)=ifixu(nod);
        u(nod)=u(nod);
    else
        ifixu(nod) = 1;
        u(nod) = 0;
    end
end
% ifixu(1) = 1;
% u(1) = 0;