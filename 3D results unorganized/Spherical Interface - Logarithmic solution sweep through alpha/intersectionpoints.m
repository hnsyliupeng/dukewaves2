function [xint,yint,zint,coeff,nb_int,edge_id] = intersectionpoints(node,x,y,z,ls,e)

xe = x(node(1:4,e));
ye = y(node(1:4,e));
ze = z(node(1:4,e));
ones = [1;1;1;1];
xmat = [ones,xe,ye,ze];
nlink = 4;
%determine which two sides have a change of sign in level-set
%Calculation of intersection points
nb_int=0;
lse = [ls(node(1,e));ls(node(2,e));ls(node(3,e));ls(node(4,e))];
coeff = xmat\lse;
for i=2:4
    for k=1:i-1
        c=lse(i)*lse(k);
        if (c<0)
            nb_int=nb_int+1;
            if(i==2)
                edge_id(nb_int)=21;
            elseif(i==3 & k==1)
                edge_id(nb_int)=31;
            elseif(i==3 & k==2)
                edge_id(nb_int)=32;
            elseif(i==4 & k==1)
                edge_id(nb_int)=41;
            elseif(i==4 & k==2)
                edge_id(nb_int)=42;
            elseif(i==4 & k==3)
                edge_id(nb_int)=43;
            end
            xint(nb_int)=xe(k)-(xe(i)-xe(k))*lse(k)/(lse(i)-lse(k));
            yint(nb_int)=ye(k)-(ye(i)-ye(k))*lse(k)/(lse(i)-lse(k));
            zint(nb_int)=ze(k)-(ze(i)-ze(k))*lse(k)/(lse(i)-lse(k));
        end
    end
end