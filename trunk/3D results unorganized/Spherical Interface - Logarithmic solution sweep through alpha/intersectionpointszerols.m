function [xint,yint,zint,coeff,nb_int,edge_id] = intersectionpointszerols(node,x,y,z,ls,e)

xe = x(node(1:4,e));
ye = y(node(1:4,e));
ze = z(node(1:4,e));
lse = ls(node(1:4,e));
ones = [1;1;1;1];
xmat = [ones,xe,ye,ze];
nlink = 4;
lstemp = zeros(nlink,1);
coeff = xmat\lse';
nb_int=0;
nb_zero = 0;
idcount = 0;
nb_pos = 0;
nb_neg = 0;
for i = 1:nlink
    if(lse(i)==0)
        nb_int = nb_int+1;
        nb_zero = nb_zero + 1;
        xint(nb_int)=xe(i); yint(nb_int)=ye(i); zint(nb_int)=ze(i);
    elseif(lse(i)>0)
        nb_pos = nb_pos + 1;
    end
end
nb_neg = nlink - nb_pos - nb_zero;
for i=2:4
    for k=1:i-1
        c=lse(i)*lse(k);
        if (c<0)
            nb_int=nb_int+1;
            xint(nb_int)=xe(k)-(xe(i)-xe(k))*lse(k)/(lse(i)-lse(k));
            yint(nb_int)=ye(k)-(ye(i)-ye(k))*lse(k)/(lse(i)-lse(k));
            zint(nb_int)=ze(k)-(ze(i)-ze(k))*lse(k)/(lse(i)-lse(k));
        end
    end
end
if (nb_int==3)
    if (nb_zero<=2)
        for i = 1:nlink
            if (lse(i)~=0)
                lstemp(i) = lse(i);
            elseif (nb_pos>nb_neg)
                lstemp(i) = lse(i)+.01;
            else
                lstemp(i) = lse(i)-.01;
            end
        end
        for i=2:4
            for k=1:i-1
                c=lstemp(i)*lstemp(k);
                if (c<0)
                    idcount=idcount+1;
                    if(i==2)
                        edge_id(idcount)=21;
                    elseif(i==3 & k==1)
                        edge_id(idcount)=31;
                    elseif(i==3 & k==2)
                        edge_id(idcount)=32;
                    elseif(i==4 & k==1)
                        edge_id(idcount)=41;
                    elseif(i==4 & k==2)
                        edge_id(idcount)=42;
                    elseif(i==4 & k==3)
                        edge_id(idcount)=43;
                    end
                end
            end
        end
    else
        count = 0;
        for i=2:4
            for k=1:i-1
                c=lse(i)*lse(k);
                if(lse(i)==0 && lse(k)==0)
                    continue;
                end
                if (c==0)
                    count=count+1;
                    if(i==2)
                        edge_id(count)=21;
                    elseif(i==3 & k==1)
                        edge_id(count)=31;
                    elseif(i==3 & k==2)
                        edge_id(count)=32;
                    elseif(i==4 & k==1)
                        edge_id(count)=41;
                    elseif(i==4 & k==2)
                        edge_id(count)=42;
                    elseif(i==4 & k==3)
                        edge_id(count)=43;
                    end
                end
            end
        end
    end
else
    edge_id = zeros(nb_int,1);
end