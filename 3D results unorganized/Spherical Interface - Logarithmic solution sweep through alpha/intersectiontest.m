function [flag1,flag2] = intersectiontest(ls,node,e)
flag1 = 0;
flag2 = 0;
for i=2:4
    schange=ls(node(i-1,e))*ls(node(i,e));
    if (schange<=0)
        flag1 = 1;
        if (schange==0)
            flag2 = 1;
        end
    end
end