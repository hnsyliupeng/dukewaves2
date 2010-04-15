function assign_parentelem_recursive(elemi,parentid)
global ELEMINFO_ARR PARENTELEM_INFO;

if (ELEMINFO_ARR(elemi).nb_subelts == 1)
    PARENTELEM_INFO(elemi) = parentid;
else
    for j=1:ELEMINFO_ARR(elemi).nb_subelts
        elemj = ELEMINFO_ARR(elemi).subelemids(j);
        assign_parentelem_recursive(elemj,parentid);
    end
end