function  [ifixu, u] = volcheck(ls,e,ifixu,u,node,volumep,tol);

lse = ls(node(1:4,e));
if(volumep<tol)
    for nlink=1:4
        if(lse(nlink)<0)
            ifixu(node(nlink,e))=1;
            u(node(nlink,e)) = 0;
        end
    end
else
    for nlink=1:4
        if(lse(nlink)<0)
            ifixu(node(nlink,e))=0;
            u(node(nlink,e)) = 0;
        end
    end
end