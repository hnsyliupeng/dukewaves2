function assign_area_to_nodes_recursive(elemi,fnodes)

global NODEINFO_ARR
global ELEMINFO_ARR CONN X Y;
global SUBELEMENT_GRAIN_MAP;

% father nodes

f1 = fnodes(1);
f2 = fnodes(2);
f3 = fnodes(3);

if (ELEMINFO_ARR(elemi).nb_subelts == 1)
    
    % Get the grain number assiciated with this (sub)element
    gnum = SUBELEMENT_GRAIN_MAP(elemi);

    if gnum == -1
        error('trying to calculate area of cut element');
    else
        %calculate area of (sub)element
        el_area = polyarea(X(CONN(1:3,elemi)),Y(CONN(1:3,elemi)));
    end

    % add the information to the nodal info
    NODEINFO_ARR(f1).areas(2,gnum) =...
        NODEINFO_ARR(f1).areas(2,gnum) + el_area;
    NODEINFO_ARR(f2).areas(2,gnum) =...
        NODEINFO_ARR(f2).areas(2,gnum) + el_area;
    NODEINFO_ARR(f3).areas(2,gnum) =...
        NODEINFO_ARR(f3).areas(2,gnum) + el_area;
    
else
    for j=1:ELEMINFO_ARR(elemi).nb_subelts
%        NODEINFO_ARR(f1).multi_grains = 1;
%        NODEINFO_ARR(f2).multi_grains = 1;
%        NODEINFO_ARR(f3).multi_grains = 1;
        elemj = ELEMINFO_ARR(elemi).subelemids(j);
        assign_area_to_nodes_recursive(elemj,fnodes)
    end
end


    