function [volumep,volumeele] = posvolume(node,e,edge_id,nb_int,x,y,z,ls,xint,yint,zint)

xe = x(node(1:4,e));
ye = y(node(1:4,e));
ze = z(node(1:4,e));
lse = ls(node(1:4,e));
volumeele = tetraunitvolume('DUMMY')*parallelipipedvolume3d(xe, ye, ze);

if (nb_int==3)
    if(edge_id(1)==21 & edge_id(2)==31 & edge_id(3)==41)
        Xt4=[xe(1),xint(1),xint(2),xint(3)]; Yt4=[ye(1),yint(1),yint(2),yint(3)]; Zt4=[ze(1),zint(1),zint(2),zint(3)];
        volumetet = tetraunitvolume('DUMMY')*parallelipipedvolume3d(Xt4, Yt4, Zt4);
        if (lse(1)>0)
            volumep = volumetet;
        elseif (lse(1)<0)
            volumep = volumeele - volumetet;
        end
    elseif(edge_id(1)==21 & edge_id(2)==32 & edge_id(3)==42)
        Xt4=[xe(2),xint(1),xint(2),xint(3)]; Yt4=[ye(2),yint(1),yint(2),yint(3)]; Zt4=[ze(2),zint(1),zint(2),zint(3)];
        volumetet = tetraunitvolume('DUMMY')*parallelipipedvolume3d(Xt4, Yt4, Zt4);
        if (lse(2)>0)
            volumep = volumetet;
        elseif (lse(2)<0)
            volumep = volumeele - volumetet;
        end
    elseif(edge_id(1)==31 & edge_id(2)==32 & edge_id(3)==43)
        Xt4=[xe(3),xint(1),xint(2),xint(3)]; Yt4=[ye(3),yint(1),yint(2),yint(3)]; Zt4=[ze(3),zint(1),zint(2),zint(3)];
        volumetet = tetraunitvolume('DUMMY')*parallelipipedvolume3d(Xt4, Yt4, Zt4);
        if (lse(3)>0)
            volumep = volumetet;
        elseif (lse(3)<0)
            volumep = volumeele - volumetet;
        end
    elseif(edge_id(1)==41 & edge_id(2)==42 & edge_id(3)==43)
        Xt4=[xe(4),xint(1),xint(2),xint(3)]; Yt4=[ye(4),yint(1),yint(2),yint(3)]; Zt4=[ze(4),zint(1),zint(2),zint(3)];
        volumetet = tetraunitvolume('DUMMY')*parallelipipedvolume3d(Xt4, Yt4, Zt4);
        if (lse(4)>0)
            volumep = volumetet;
        elseif (lse(4)<0)
            volumep = volumeele - volumetet;
        end
    end

    %% Four intersection points
elseif (nb_int==4)
    if(edge_id(1)==21 & edge_id(2)==32 & edge_id(3)==41 & edge_id(4)==43)
        Xt1=[xe(1),xint(3),xint(2),xe(3)]; Yt1=[ye(1),yint(3),yint(2),ye(3)]; Zt1=[ze(1),zint(3),zint(2),ze(3)];
        volumet1 = tetraunitvolume('DUMMY')*parallelipipedvolume3d(Xt1, Yt1, Zt1);

        Xt2=[xint(4),xint(3),xe(3),xint(2)]; Yt2=[yint(4),yint(3),ye(3),yint(2)]; Zt2=[zint(4),zint(3),ze(3),zint(2)];
        volumet2 = tetraunitvolume('DUMMY')*parallelipipedvolume3d(Xt2, Yt2, Zt2);

        Xt3=[xint(3),xe(1),xint(2),xint(1)]; Yt3=[yint(3),ye(1),yint(2),yint(1)]; Zt3=[zint(3),ze(1),zint(2),zint(1)];
        volumet3 = tetraunitvolume('DUMMY')*parallelipipedvolume3d(Xt3, Yt3, Zt3);

        volumewedge = volumet1 + volumet2 + volumet3;
        
        if(lse(1) && lse(3)>0)
            volumep = volumewedge;
        else
            volumep = volumeele - volumewedge;
        end

    elseif(edge_id(1)==21 & edge_id(2)==31 & edge_id(3)==42 & edge_id(4)==43)
        Xt1=[xe(4),xint(1),xint(4),xe(1)]; Yt1=[ye(4),yint(1),yint(4),ye(1)]; Zt1=[ze(4),zint(1),zint(4),ze(1)];
        volumet1 = tetraunitvolume('DUMMY')*parallelipipedvolume3d(Xt1, Yt1, Zt1);

        Xt2=[xint(3),xint(4),xe(4),xint(1)]; Yt2=[yint(3),yint(4),ye(4),yint(1)]; Zt2=[zint(3),zint(4),ze(4),zint(1)];
        volumet2 = tetraunitvolume('DUMMY')*parallelipipedvolume3d(Xt2, Yt2, Zt2);

        Xt3=[xint(4),xe(1),xint(1),xint(2)]; Yt3=[yint(4),ye(1),yint(1),yint(2)]; Zt3=[zint(4),ze(1),zint(1),zint(2)];
        volumet3 = tetraunitvolume('DUMMY')*parallelipipedvolume3d(Xt3, Yt3, Zt3);

        volumewedge = volumet1 + volumet2 + volumet3;
        
        if(lse(1) && lse(4)>0)
            volumep = volumewedge;
        else
            volumep = volumeele - volumewedge;
        end

    elseif(edge_id(1)==31 & edge_id(2)==32 & edge_id(3)==41 & edge_id(4)==42)
        Xt1=[xe(1),xint(1),xint(4),xe(2)]; Yt1=[ye(1),yint(1),yint(4),ye(2)]; Zt1=[ze(1),zint(1),zint(4),ze(2)];
        volumet1 = tetraunitvolume('DUMMY')*parallelipipedvolume3d(Xt1, Yt1, Zt1);

        Xt2=[xint(2),xint(4),xe(2),xint(1)]; Yt2=[yint(2),yint(4),ye(2),yint(1)]; Zt2=[zint(2),zint(4),ze(2),zint(1)];
        volumet2 = tetraunitvolume('DUMMY')*parallelipipedvolume3d(Xt2, Yt2, Zt2);

        Xt3=[xint(4),xe(1),xint(1),xint(3)]; Yt3=[yint(4),ye(1),yint(1),yint(3)]; Zt3=[zint(4),ze(1),zint(1),zint(3)];
        volumet3 = tetraunitvolume('DUMMY')*parallelipipedvolume3d(Xt3, Yt3, Zt3);

        volumewedge = volumet1 + volumet2 + volumet3;

        if(lse(1) && lse(2)>0)
            volumep = volumewedge;
        else
            volumep = volumeele - volumewedge;
        end
       
    end
end