function [ L2_error, L2_exact, domvol] = L2errorcutelem(node,e,B,x,y,z,xint,yint,zint,disp,numnod,ls,ifixu,nb_int,edge_id,volumep);

xe = x(node(1:4,e));
ye = y(node(1:4,e));
ze = z(node(1:4,e));
ue = disp(node(1:4,e));
lse = ls(node(1:4,e));
Amat = [ 1 1 1 1; xe(1) xe(2) xe(3) xe(4); ye(1) ye(2) ye(3) ye(4); ze(1) ze(2) ze(3) ze(4)];
volume_sub = volumep;
wedgeflag = 0;
if (nb_int==3)
    if(edge_id(1)==21 & edge_id(2)==31 & edge_id(3)==41)
        %%%%%% Subelement is a tetrahedron %%%%%%%%
        if (lse(1)>0)
            Xt4=[xe(1),xint(1),xint(2),xint(3)]; Yt4=[ye(1),yint(1),yint(2),yint(3)]; Zt4=[ze(1),zint(1),zint(2),zint(3)];
            %%%%%% Subelement is a wedge %%%%%%%%%%%%%
        elseif (lse(1)<0)
            wedgeflag = 1;
        end
    elseif(edge_id(1)==21 & edge_id(2)==32 & edge_id(3)==42)
        %%%%%%%% Subelement is a terahedron %%%%%%%%%%%%%%%
        if (lse(2)>0)
            Xt4=[xe(2),xint(1),xint(2),xint(3)]; Yt4=[ye(2),yint(1),yint(2),yint(3)]; Zt4=[ze(2),zint(1),zint(2),zint(3)];
            %%%%%%%%% Subelement is a wedge %%%%%%%%%%%%%%%%%%
        elseif (lse(2)<0)
            wedgeflag = 1;
        end
    elseif(edge_id(1)==31 & edge_id(2)==32 & edge_id(3)==43)
        %%%%%%%%%% Subelement is a tetrahedron %%%%%%%%%%%%%
        if (lse(3)>0)
            Xt4=[xe(3),xint(1),xint(2),xint(3)]; Yt4=[ye(3),yint(1),yint(2),yint(3)]; Zt4=[ze(3),zint(1),zint(2),zint(3)];
            %%%%%%%%%% Subelement is a wedge %%%%%%%%%%%%%%%%%%
        elseif (lse(3)<0)
            wedgeflag = 1;
        end
    elseif(edge_id(1)==41 & edge_id(2)==42 & edge_id(3)==43)
        %%%%%%%%%% Subelement is a tetrahedron %%%%%%%%%%%%%%%%%
        if (lse(4)>0)
            Xt4=[xe(4),xint(1),xint(2),xint(3)]; Yt4=[ye(4),yint(1),yint(2),yint(3)]; Zt4=[ze(4),zint(1),zint(2),zint(3)];
            %%%%%%%%%% Subelement is a wedge %%%%%%%%%%%%%%%%%%%%%%%
        elseif (lse(4)<0)
            wedgeflag = 1;
        end
    end
    %%%%%%% Tetrahedral sub-element %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(wedgeflag==0)
        xe_sub = Xt4'; ye_sub = Yt4'; ze_sub = Zt4';
        [L2_error,L2_exact, domvol] = L2errortetsub(xe_sub,ye_sub,ze_sub,volume_sub,Amat,ue);
    else
        %%%%%%%%%%%%% Wedge Shaped sub-element %%%%%%%%%%%%%%%%%%%%%%%%
        if(edge_id(1)==21 & edge_id(2)==31 & edge_id(3)==41)
            xe_sub = [xe(2);xe(3);xe(4);xint(1);xint(2);xint(3)];
            ye_sub = [ye(2);ye(3);ye(4);yint(1);yint(2);yint(3)];
            ze_sub = [ze(2);ze(3);ze(4);zint(1);zint(2);zint(3)];
            vector_41 = [xe(2)-xe(1),ye(2)-ye(1),ze(2)-ze(1)];
            vector_52 = [xe(3)-xe(1),ye(3)-ye(1),ze(3)-ze(1)];
            [xe_sub,ye_sub,ze_sub] = orderwedgenodes(xe_sub,ye_sub,ze_sub,xe,ye,ze,vector_41,vector_52);
        elseif(edge_id(1)==21 & edge_id(2)==32 & edge_id(3)==42)
            xe_sub = [xe(1);xe(3);xe(4);xint(1);xint(2);xint(3)];
            ye_sub = [ye(1);ye(3);ye(4);yint(1);yint(2);yint(3)];
            ze_sub = [ze(1);ze(3);ze(4);zint(1);zint(2);zint(3)];
            vector_41 = [xe(2)-xe(1),ye(2)-ye(1),ze(2)-ze(1)];
            vector_52 = [xe(2)-xe(3),ye(2)-ye(3),ze(2)-ze(3)];
            [xe_sub,ye_sub,ze_sub] = orderwedgenodes(xe_sub,ye_sub,ze_sub,xe,ye,ze,vector_41,vector_52);
        elseif(edge_id(1)==31 & edge_id(2)==32 & edge_id(3)==43)
            xe_sub = [xe(1);xe(2);xe(4);xint(1);xint(2);xint(3)];
            ye_sub = [ye(1);ye(2);ye(4);yint(1);yint(2);yint(3)];
            ze_sub = [ze(1);ze(2);ze(4);zint(1);zint(2);zint(3)];
            vector_41 = [xe(3)-xe(1),ye(3)-ye(1),ze(3)-ze(1)];
            vector_52 = [xe(3)-xe(2),ye(3)-ye(2),ze(3)-ze(2)];
            [xe_sub,ye_sub,ze_sub] = orderwedgenodes(xe_sub,ye_sub,ze_sub,xe,ye,ze,vector_41,vector_52);
        elseif(edge_id(1)==41 & edge_id(2)==42 & edge_id(3)==43)
            xe_sub = [xe(1);xe(2);xe(3);xint(1);xint(2);xint(3)];
            ye_sub = [ye(1);ye(2);ye(3);yint(1);yint(2);yint(3)];
            ze_sub = [ze(1);ze(2);ze(3);zint(1);zint(2);zint(3)];
            vector_41 = [xe(4)-xe(1),ye(4)-ye(1),ze(4)-ze(1)];
            vector_52 = [xe(4)-xe(2),ye(4)-ye(2),ze(4)-ze(2)];
            [xe_sub,ye_sub,ze_sub] = orderwedgenodes(xe_sub,ye_sub,ze_sub,xe,ye,ze,vector_41,vector_52);
        end
%         if(e==64)
%             xe_sub = [0;1;0;0;0.35;0];
%             ye_sub = [0;1;1;0;0.35;0.35];
%             ze_sub = [11/23;11/23;11/23;.45;.45;.45];
%         elseif(e==63)
%             xe_sub = [1;1;0;0.35;0.35;0];
%             ye_sub = [0;1;0;0;0.35;0];
%             ze_sub = [11/23;11/23;11/23;.45;.45;.45];
%         end
        [L2_error, L2_exact, domvol] = L2errorwedgesub(xe_sub,ye_sub,ze_sub,volume_sub,Amat,ue);
        %%%%% Wedge Routine Verification %%%%%%%%%%%%%%%%%%
        [wedgevolume] = wedgesubveri(xe_sub,ye_sub,ze_sub,Amat);
    end
elseif(nb_int==4)
    if(edge_id(1)==21 & edge_id(2)==32 & edge_id(3)==41 & edge_id(4)==43)
        if(lse(1) && lse(3)>0)
            %%%%%%% Top-wedge is in the domain %%%%%%%%%%%%%%%%%%
            xe_sub = [xint(1);xint(2);xe(1);xint(4);xint(3);xe(3)];
            ye_sub = [yint(1);yint(2);ye(1);yint(4);yint(3);ye(3)];
            ze_sub = [zint(1);zint(2);ze(1);zint(4);zint(3);ze(3)];
            vector_31 = [xe(2)-xe(1),ye(2)-ye(1),ze(2)-ze(1)];  %% vector_sub_31 should lie on vector_31
            vector_32 = [xe(4)-xe(1),ye(4)-ye(1),ze(4)-ze(1)];  %% vector_sub_32 should lie on vector_32
            [xe_sub,ye_sub,ze_sub] = orderfourintwedgenodes(xe_sub,ye_sub,ze_sub,vector_31,vector_32);
        else
            %%%%%%%%% Bottom wedge is in the domain %%%%%%%%%%%%%%%
            xe_sub = [xint(1);xint(4);xe(2);xint(2);xint(3);xe(4)];
            ye_sub = [yint(1);yint(4);ye(2);yint(2);yint(3);ye(4)];
            ze_sub = [zint(1);zint(4);ze(2);zint(2);zint(3);ze(4)];
            vector_31 = [xe(2)-xe(1),ye(2)-ye(1),ze(2)-ze(1)];
            vector_32 = [xe(3)-xe(2),ye(3)-ye(2),ze(3)-ze(2)];
            [xe_sub,ye_sub,ze_sub] = orderfourintwedgenodes(xe_sub,ye_sub,ze_sub,vector_31,vector_32);
        end

    elseif(edge_id(1)==21 & edge_id(2)==31 & edge_id(3)==42 & edge_id(4)==43)
        if(lse(1) && lse(4)>0)
            %%%%%%% Top-wedge is in the domain %%%%%%%%%%%%%%%%%%
            xe_sub = [xint(1);xint(2);xe(1);xint(4);xint(3);xe(4)];
            ye_sub = [yint(1);yint(2);ye(1);yint(4);yint(3);ye(4)];
            ze_sub = [zint(1);zint(2);ze(1);zint(4);zint(3);ze(4)];
            vector_31 = [xe(2)-xe(1),ye(2)-ye(1),ze(2)-ze(1)];
            vector_32 = [xe(3)-xe(1),ye(3)-ye(1),ze(3)-ze(1)];
            [xe_sub,ye_sub,ze_sub] = orderfourintwedgenodes(xe_sub,ye_sub,ze_sub,vector_31,vector_32);
        else
            %%%%%%%%% Bottom wedge is in the domain %%%%%%%%%%%%%%%
            xe_sub = [xint(1);xint(4);xe(2);xint(2);xint(3);xe(3)];
            ye_sub = [yint(1);yint(4);ye(2);yint(2);yint(3);ye(3)];
            ze_sub = [zint(1);zint(4);ze(2);zint(2);zint(3);ze(3)];
            vector_31 = [xe(2)-xe(1),ye(2)-ye(1),ze(2)-ze(1)];
            vector_32 = [xe(2)-xe(4),ye(2)-ye(4),ze(2)-ze(4)];
            [xe_sub,ye_sub,ze_sub] = orderfourintwedgenodes(xe_sub,ye_sub,ze_sub,vector_31,vector_32);
        end

    elseif(edge_id(1)==31 & edge_id(2)==32 & edge_id(3)==41 & edge_id(4)==42)
        if(lse(1) && lse(2)>0)
            %%%%%%% Top-wedge is in the domain %%%%%%%%%%%%%%%%%%
            xe_sub = [xint(1);xint(2);xe(1);xint(4);xint(3);xe(2)];
            ye_sub = [yint(1);yint(2);ye(1);yint(4);yint(3);ye(2)];
            ze_sub = [zint(1);zint(2);ze(1);zint(4);zint(3);ze(2)];
            vector_31 = [xe(3)-xe(1),ye(3)-ye(1),ze(3)-ze(1)];
            vector_32 = [xe(4)-xe(1),ye(4)-ye(1),ze(4)-ze(1)];
            [xe_sub,ye_sub,ze_sub] = orderfourintwedgenodes(xe_sub,ye_sub,ze_sub,vector_31,vector_32);
        else
            %%%%%%%%% Bottom wedge is in the domain %%%%%%%%%%%%%%%
            xe_sub = [xint(1);xint(4);xe(3);xint(2);xint(3);xe(4)];
            ye_sub = [yint(1);yint(4);ye(3);yint(2);yint(3);ye(4)];
            ze_sub = [zint(1);zint(4);ze(3);zint(2);zint(3);ze(4)];
            vector_31 = [xe(3)-xe(1),ye(3)-ye(1),ze(3)-ze(1)];
            vector_32 = [xe(3)-xe(2),ye(3)-ye(2),ze(3)-ze(2)];
            [xe_sub,ye_sub,ze_sub] = orderfourintwedgenodes(xe_sub,ye_sub,ze_sub,vector_31,vector_32);
        end
    end
    %%%%%%%%%%% Wedge Routine Verification %%%%%%%%%%%%%%
    [wedgevolume] = wedgesubveri(xe_sub,ye_sub,ze_sub,Amat);
    [L2_error, L2_exact, domvol] = L2errorwedgesub(xe_sub,ye_sub,ze_sub,volume_sub,Amat,ue);
end