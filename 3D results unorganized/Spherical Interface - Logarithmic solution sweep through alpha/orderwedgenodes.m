function [xe_sub,ye_sub,ze_sub] = orderwedgenodes(xe_sub,ye_sub,ze_sub,xe,ye,ze,vector_41,vector_52);

flag=0;
count=0;
tol = 1.0e-06;

while(flag==0)
    
    modvector_41 = sqrt(vector_41(1)*vector_41(1) + vector_41(2)*vector_41(2) + vector_41(3)*vector_41(3));
    unitvector_41 = vector_41/modvector_41;
    vector_sub_41 = [xe_sub(4)-xe_sub(1),ye_sub(4)-ye_sub(1),ze_sub(4)-ze_sub(1)];
    modvector_sub_41 = sqrt(vector_sub_41(1)*vector_sub_41(1) + vector_sub_41(2)*vector_sub_41(2) + vector_sub_41(3)*vector_sub_41(3));
    unitvector_sub_41 = vector_sub_41/modvector_sub_41;
    directioncheck_1 = dot(unitvector_41,unitvector_sub_41);
    checkvar1 = abs(directioncheck_1-1);
    checkvar2 = abs(directioncheck_1+1);

    if(checkvar1<tol || checkvar2<tol)
        flag = 1;
    else
        count = count + 1;
        if(count==1)
            xtemp = xe_sub(4); ytemp = ye_sub(4); ztemp = ze_sub(4);
            xe_sub(4) = xe_sub(5); ye_sub(4) = ye_sub(5); ze_sub(4) = ze_sub(5);
            xe_sub(5) = xtemp; ye_sub(5) = ytemp; ze_sub(5) = ztemp;
        elseif(count==2)
            xtemp = xe_sub(4); ytemp = ye_sub(4); ztemp = ze_sub(4);
            xe_sub(4) = xe_sub(6); ye_sub(4) = ye_sub(5); ze_sub(4) = ze_sub(5);
            xe_sub(5) = xtemp; ye_sub(5) = ytemp; ze_sub(5) = ztemp;
        end
    end
end
    

modvector_52 = sqrt(vector_52(1)*vector_52(1) + vector_52(2)*vector_52(2) + vector_52(3)*vector_52(3));
unitvector_52 = vector_52/modvector_52;
vector_sub_52 = [xe_sub(5)-xe_sub(2),ye_sub(5)-ye_sub(2),ze_sub(5)-ze_sub(2)];
modvector_sub_52 = sqrt(vector_sub_52(1)*vector_sub_52(1) + vector_sub_52(2)*vector_sub_52(2) + vector_sub_52(3)*vector_sub_52(3));
unitvector_sub_52 = vector_sub_52/modvector_sub_52;
directioncheck_2 = dot(unitvector_52,unitvector_sub_52);

checkvar3 = abs(directioncheck_2-1);
checkvar4 = abs(directioncheck_2+1);
flag2=0;
if(checkvar3<tol || checkvar4<tol)
    flag2 = 1;
end
if(flag2==0)
    xtemp = xe_sub(5); ytemp = ye_sub(5); ztemp = ze_sub(5);
    xe_sub(5) = xe_sub(6); ye_sub(5) = ye_sub(6); ze_sub(5) = ze_sub(6);
    xe_sub(6) = xtemp; ye_sub(6) = ytemp; ze_sub(6) = ztemp;
end