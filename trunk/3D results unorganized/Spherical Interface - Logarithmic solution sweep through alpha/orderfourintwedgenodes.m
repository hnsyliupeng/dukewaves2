function [xe_sub,ye_sub,ze_sub] = orderfourintwedgenodes(xe_sub,ye_sub,ze_sub,vector_31,vector_32);

flag=0;
count=0;
tol = 1.0e-06;

while(flag==0)

    modvector_31 = sqrt(vector_31(1)*vector_31(1) + vector_31(2)*vector_31(2) + vector_31(3)*vector_31(3));
    unitvector_31 = vector_31/modvector_31;
    vector_sub_31 = [xe_sub(3)-xe_sub(1),ye_sub(3)-ye_sub(1),ze_sub(3)-ze_sub(1)];
    modvector_sub_31 = sqrt(vector_sub_31(1)*vector_sub_31(1) + vector_sub_31(2)*vector_sub_31(2) + vector_sub_31(3)*vector_sub_31(3));
    unitvector_sub_31 = vector_sub_31/modvector_sub_31;
    directioncheck_1 = dot(unitvector_31,unitvector_sub_31);
    checkvar1 = abs(directioncheck_1 - 1);
    checkvar2 = abs(directioncheck_1 + 1);

    if(checkvar1<tol || checkvar2<tol)
        flag = 1;
    else
        for count=1:3
            if(count==1)
                xtemp = xe_sub(4); ytemp = ye_sub(4); ztemp = ze_sub(4);
                xe_sub(4) = xe_sub(1); ye_sub(4) = ye_sub(1); ze_sub(4) = ze_sub(1);
                xe_sub(1) = xtemp; ye_sub(1) = ytemp; ze_sub(1) = ztemp;
            elseif(count==2)
                xtemp = xe_sub(4); ytemp = ye_sub(4); ztemp = ze_sub(4);
                xe_sub(4) = xe_sub(2); ye_sub(4) = ye_sub(2); ze_sub(4) = ze_sub(2);
                xe_sub(2) = xtemp; ye_sub(2) = ytemp; ze_sub(2) = ztemp;
            elseif(count==3)
                xtemp = xe_sub(4); ytemp = ye_sub(4); ztemp = ze_sub(4);
                xe_sub(4) = xe_sub(5); ye_sub(4) = ye_sub(5); ze_sub(4) = ze_sub(5);
                xe_sub(5) = xtemp; ye_sub(5) = ytemp; ze_sub(5) = ztemp;
            end
        end
    end
    modvector_32 = sqrt(vector_32(1)*vector_32(1) + vector_32(2)*vector_32(2) + vector_32(3)*vector_32(3));
    unitvector_32 = vector_32/modvector_32;
    vector_sub_32 = [xe_sub(3)-xe_sub(2),ye_sub(3)-ye_sub(2),ze_sub(3)-ze_sub(2)];
    modvector_sub_32 = sqrt(vector_sub_32(1)*vector_sub_32(1) + vector_sub_32(2)*vector_sub_32(2) + vector_sub_32(3)*vector_sub_32(3));
    unitvector_sub_32 = vector_sub_32/modvector_sub_32;
    directioncheck_2 = dot(unitvector_32,unitvector_sub_32);
    checkvar3 = abs(directioncheck_2-1);
    checkvar4 = abs(directioncheck_2+1);

    if(checkvar3<tol || checkvar4<tol )
        flag = 1;
    else
        flag = 0;
        xtemp = xe_sub(4); ytemp = ye_sub(4); ztemp = ze_sub(4);
        xe_sub(4) = xe_sub(2); ye_sub(4) = ye_sub(2); ze_sub(4) = ze_sub(2);
        xe_sub(2) = xtemp; ye_sub(2) = ytemp; ze_sub(2) = ztemp;
    end
end

xe_sub = xe_sub; ye_sub = ye_sub; ze_sub = ze_sub;