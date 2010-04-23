function assign_grains_to_elements_recursive(p, tri, elemi)
global ELEMINFO_ARR CONN X Y;
global SUBELEMENT_GRAIN_MAP;

if (ELEMINFO_ARR(elemi).nb_subelts == 1)
    %find the center of the element
    xcen = (X(CONN(1,elemi)) + X(CONN(2,elemi)) + X(CONN(3,elemi)) )/3.0;
    ycen = (Y(CONN(1,elemi)) + Y(CONN(2,elemi)) + Y(CONN(3,elemi)) )/3.0;
    k = dsearch(p(:,1), p(:,2), tri, xcen, ycen);
    SUBELEMENT_GRAIN_MAP(elemi) = k;
%     figure(3)
%    patch(X(CONN(1:3,elemi)),Y(CONN(1:3,elemi)),k)
%    drawnow
else
    for j=1:ELEMINFO_ARR(elemi).nb_subelts
        elemj = ELEMINFO_ARR(elemi).subelemids(j);
        assign_grains_to_elements_recursive(p, tri, elemj);
    end
end