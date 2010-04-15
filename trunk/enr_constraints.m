function [ke_con,id_node,id_lag,rhs] =...
    enr_constraints(node,x,y,i,id_eqns,id_dof,enr_surf,temp,ubar,dispbc,nodegrainmap)

    % INITIALIZE
    
    surf = i;
    
    xs = [];
    ys = [];
    
    xsi = [];
    eta = [];
    
    N = zeros(2,1);
    
    g = zeros(1,8);

    % Get node numbers
    n1 = enr_surf(i).nodes(1);
    n2 = enr_surf(i).nodes(2);
    
    % Get coordinates of nodes
    xs = x([n1 n2]);
    ys = y([n1 n2]);
    
    % Get segment jacobian
    p1 = enr_surf(i).coords(1,:);
    p2 = enr_surf(i).coords(2,:);
    
    % jacobian of segment to global
    he = sqrt((p1(1)-p2(1))^2 + (p1(2)-p2(2))^2);
    seg_jcob = he/2;    

    % Gauss points on segments
    gauss = [-sqrt(3)/3 sqrt(3)/3];
    weights = [1 1];
        
    eta = [-1 1];
    
    xsi = enr_surf(i).xsi;

    % loop over Gauss points
    
    for g = 1:2
        
        xn = 0.5*(1-gauss(g))*xsi(1)+0.5*(1+gauss(g))*xsi(2);
    
        for b = 1:2    % Evaluate 2 shape functions        
            N(b) = N(b) + 0.5*(1 + eta(b)*xn)*seg_jcob*weights(g);        
        end
    end

    ke_con = zeros(2,8);
    ke_con(1,1) = N(1);
    ke_con(2,2) = N(1);
    ke_con(1,5) = N(2);
    ke_con(2,6) = N(2);
    if enr_surf(i).grain == id_dof(enr_surf(i).nodes(1),3)
        ke_con(1,3) = N(1);
        ke_con(2,4) = N(1);
    end
    if enr_surf(i).grain == id_dof(enr_surf(i).nodes(2),3)
        ke_con(1,7) = N(2);
        ke_con(2,8) = N(2);
    end
    
    id_lag = [temp+2*surf-1, temp+2*surf];
    id_node(1) = id_eqns(n1,1);
    id_node(2) = id_eqns(n1,2);
    id_node(3) = id_eqns(n1,3);
    id_node(4) = id_eqns(n1,4);    
    id_node(5) = id_eqns(n2,1);
    id_node(6) = id_eqns(n2,2);
    id_node(7) = id_eqns(n2,3);
    id_node(8) = id_eqns(n2,4); 
    
    % Assemble right hand side
    
    g(1) = ubar(1,n1);
    g(2) = ubar(2,n1);
    g(3) = 0;
    g(4) = 0;
    g(5) = ubar(1,n2);
    g(6) = ubar(2,n2);
    g(7) = 0;
    g(8) = 0;
    rhs = ke_con*g';
