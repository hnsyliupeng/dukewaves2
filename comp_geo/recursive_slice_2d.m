function recursive_slice_2d(elemi, cuteleminfo)

global ELEMINFO_ARR;

%root of recursive algorithm
if (ELEMINFO_ARR(elemi).nb_subelts == 1)
    [nb_sub_elts, sub_elt_ids] = slice_2d_element_by_segment(elemi, cuteleminfo);
    eleminfo = struct('nb_subelts',nb_sub_elts,'subelemids',sub_elt_ids); 
    ELEMINFO_ARR(elemi) = eleminfo;
    if (nb_sub_elts > 1)
       for j=1:nb_sub_elts
           eleminfo = struct('nb_subelts',1,'subelemids',-1);
           ELEMINFO_ARR(sub_elt_ids(j)) = eleminfo;
       end
    end
else
    %elemi
    %ELEMINFO_ARR(elemi).nb_subelts;
    for j=1:ELEMINFO_ARR(elemi).nb_subelts
        elemj = ELEMINFO_ARR(elemi).subelemids(j);
        recursive_slice_2d(elemj, cuteleminfo);
    end
end