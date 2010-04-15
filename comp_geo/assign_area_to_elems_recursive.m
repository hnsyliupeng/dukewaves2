function assign_area_to_elems_recursive(elemi,father)

global ELEMINFO_ARR CONN X Y;
global SUBELEMENT_GRAIN_MAP;
global FATHERS_ARR

if (ELEMINFO_ARR(elemi).nb_subelts == 1)
    
    % Get the grain number assiciated with this (sub)element
    gnum = SUBELEMENT_GRAIN_MAP(elemi);

    if gnum == -1
        error('trying to calculate area of cut element');
    else
        %calculate area of (sub)element
        el_area = polyarea(X(CONN(1:3,elemi)),Y(CONN(1:3,elemi)));
    end

    % add the area to the father element info
    FATHERS_ARR(father).areas(2,gnum) =...
        FATHERS_ARR(father).areas(2,gnum) + el_area;
    
else
    for j=1:ELEMINFO_ARR(elemi).nb_subelts
%        NODEINFO_ARR(f1).multi_grains = 1;
%        NODEINFO_ARR(f2).multi_grains = 1;
%        NODEINFO_ARR(f3).multi_grains = 1;
        elemj = ELEMINFO_ARR(elemi).subelemids(j);
        assign_area_to_elems_recursive(elemj,father)
    end
end


    