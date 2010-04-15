function [id,ke] = elemstiff_class_recursive(node,x,y,elemi,parentid,id_dof,id_eqns)
global ELEMINFO_ARR
global NODAL_ENRICH
global PARENTELEM_INFO

    % Determine here the dimension of the stiffness
    dof = zeros(1,3);

    % Determine how many degrees of freedom each node has
    for i = 1:3
        dof(i) = 2 + 2*(NODAL_ENRICH(node(i,parentid)).cnt - 1);
    end                                                             
                                                                   
    tot_dof = dof(1) + dof(2) + dof(3);                                            
    ke = zeros(tot_dof,tot_dof);

    id = [id_eqns(node(1,parentid),1:2) id_eqns(node(2,parentid),1:2) id_eqns(node(3,parentid),1:2)...
        id_eqns(node(1,parentid),3:6) id_eqns(node(2,parentid),3:6) id_eqns(node(3,parentid),3:6)];
    eliminate = find(id == 0);
    for i = size(eliminate,2):-1:1
        id(eliminate(i)) = [];
    end

if (ELEMINFO_ARR(elemi).nb_subelts == 1)
    parentid = PARENTELEM_INFO(elemi);
%    elemi
    [ke] = elemstiff_class_subelement(node,x,y,elemi,parentid,id_dof);
else
    for j=1:ELEMINFO_ARR(elemi).nb_subelts
        elemj = ELEMINFO_ARR(elemi).subelemids(j);
        [id,ke_sub] = elemstiff_class_recursive(node,x,y,elemj,parentid,id_dof,id_eqns);
        
        %add the contributions
        ke = ke + ke_sub;
    end
end
