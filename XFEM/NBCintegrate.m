% NBCintegrate.m
%
% CALL: NBCintegrate(BOUNDARY,FORCE,DOFs1,DOFs2)
% 
% Integrate over a single, uncut element edge to get the nodal forces, that
% have to be assembled into 'big_force'
%
% Input parameters
%   BOUNDARY          structure, that describes the boundary edge
%                     (integration domain)
%   FORCE             stucture, that describes the traction, that is given
%                     in this element edge
%   DOFs1             global degrees of freedom for node 1 of element edge
%   DOFs2             global degrees of freedom for node 2 of element edge
%
% Returned variables
%   force_values      element load vector (nodal forces)
%   force_id          ID array for assembly

% Author: Matthias Mayr (05/2010)

function [force_values,force_id] = NBCintegrate(BOUNDARY,FORCE,DOFs1,DOFs2)

if strcmpi(FORCE.shape,'constant')
    % compute length of Neumann boundary
    len = norm(FORCE.coords(1,:) - FORCE.coords(2,:));
    
    % constant force --> apply half of 'FORCE.values' on each node
    force_values = 0.5 * [FORCE.values(1) FORCE.values(2) ...
        FORCE.values(1) FORCE.values(2)];% ./ len;
    force_id = [DOFs1(1) DOFs1(2) DOFs2(1) DOFs2(2)];
elseif strcmp(FORCE.shape,'linear')
    % Now, use 2-point-gauss-integration to evaluate the integral over the
    % Neumman boundary: \int _{\Gamma _h} w _i h _i d \Gamma
    % with  w   test function
    %       h   traction vector
    
    % define some values for 2-point-gauss integration
    gauss = [-sqrt(3)/3 sqrt(3)/3];
    weights = [1.0 1.0];
    
    % compute jacobian
    jac = norm(BOUNDARY.coords(1,:) - BOUNDARY.coords(2,:))/2;
    
    % construct function, that describes forces: f(x) = m*x + t
    fx1 = FORCE.values(1,1);    % value of x-force on one end
    fy1 = FORCE.values(1,2);    % value of y-force on one end
    fx2 = FORCE.values(2,1);    % value of x-force on other end
    fy2 = FORCE.values(2,2);    % value of y-force on other end
    
    % compute length of Neumann boundary
    len = norm(FORCE.coords(1,:) - FORCE.coords(2,:));
    
    % slope and axis intersept of linear traction function (origin chosen
    % at one end
    mx = (fx2 - fx1) / len;
    tx = fx1;
    my = (fy2 - fy1) / len;
    ty = fy1;
    
    % define function for linear traction
    traction_x = @(x) mx*(x+len/2) + tx;
    traction_y = @(x) my*x + ty;
        
    % get coords of element nodes
    p1 = BOUNDARY.coords(1,:);
    p2 = BOUNDARY.coords(2,:);
    
    % Get real coordinates of gauss points
    xgp1 = 0.5*(1-gauss(1))*p1(1)+0.5*(1+gauss(1))*p2(1);
    ygp1 = 0.5*(1-gauss(1))*p1(2)+0.5*(1+gauss(1))*p2(2);
    xgp2 = 0.5*(1-gauss(2))*p1(1)+0.5*(1+gauss(2))*p2(1);
    ygp2 = 0.5*(1-gauss(2))*p1(2)+0.5*(1+gauss(2))*p2(2);
    
    % evaluate traction vector 'h' at gauss points
    hgp1 = [feval(traction_x,ygp1);feval(traction_y,xgp1)]; 
    hgp2 = [feval(traction_x,ygp2);feval(traction_y,xgp2)]; 
    
    % Due to linear triangular elements, every shape function along an
    % element edge is linear: N_1 = 0.5(1 - xsi), N_2 = 0.5(1 + xsi)
    % After evaluation at the gauss points 'gauss', they are the same for
    % every edge, if integrating along an element edge (Neumann boundary).
    % So, they are hard-coded here.
    N_11 = 0.5*(1-gauss(1));    % shape function 1 at gauss point 1
    N_12 = 0.5*(1+gauss(1));    % shape function 2 at gauss point 1
    N_21 = 0.5*(1-gauss(2));    % shape function 1 at gauss point 2
    N_22 = 0.5*(1+gauss(2));    % shape function 2 at gauss point 2
    
    % build shape function matrices 'Ngp1' and 'Ngp2' at gauss points 1 & 2
    Ngp1 = [N_11 0 N_12 0;
            0 N_11 0 N_12];
    Ngp2 = [N_21 0 N_22 0;
        0 N_21 0 N_22];
    
    force_values = Ngp1' * hgp1 * jac * weights(1) + ... 
        Ngp2' * hgp2 * jac * weights(2);
    force_id = [DOFs1(1) DOFs1(2) DOFs2(1) DOFs2(2)];
    
elseif strcmp(FORCE.shape,'parabolic')
    % Now, use 3-point-gauss-integration to evaluate the integral over the
    % Neumman boundary: \int _{\Gamma _h} w _i h _i d \Gamma
    % with  w   test function
    %       h   traction vector
    %
    % The parabolic traction distribution is computed via 2 given values at
    % the one end and in the middle of the parabel. Lagrange polynomials
    % are used to get a function.
    
    % define some values for 2-point-gauss integration
    gauss = [-sqrt(3/5) 0 sqrt(3/5)];
    weights = [5/9 8/9 5/9];
    
    % compute jacobian
    jac = norm(BOUNDARY.coords(1,:) - BOUNDARY.coords(2,:))/2;
    
    
    % comput length of Neumann boundary
    len = norm(FORCE.coords(1,:) - FORCE.coords(2,:));
    
    % construct function, that describes Lagrange polynomial for parabolic
    % force distribution
    fx1 = FORCE.values(1,1);    % value of x-force on one end
    fy1 = FORCE.values(1,2);    % value of y-force on one end
    fx2 = FORCE.values(2,1);    % value of x-force on other end
    fy2 = FORCE.values(2,2);    % value of y-force on other end
    
    % define function for linear traction (via Lagrange polynomials)
    traction_x = @(x) (x-0.5*len)*(x-len)/(-0.5*len)*(-len) * fx1 + ...
        x*(x-len)/(0.5*len)*(-0.5*len) * fx2 + ...
        x*(x-0.5*len)/(len)*(0.5*len) * fx1;
    traction_y = @(y) (y-0.5*len)*(y-len)/(-0.5*len)*(-len) * fy1 + ...
        y*(y-len)/(0.5*len)*(-0.5*len) * fy2 + ...
        y*(y-0.5*len)/(len)*(0.5*len) * fy1;
    
    % get coords of element nodes
    p1 = BOUNDARY.coords(1,:);
    p2 = BOUNDARY.coords(2,:);
    
    % Get real coordinates of gauss points
    xgp1 = 0.5*(1-gauss(1))*p1(1)+0.5*(1+gauss(1))*p2(1);
    ygp1 = 0.5*(1-gauss(1))*p1(2)+0.5*(1+gauss(1))*p2(2);
    xgp2 = 0.5*(1-gauss(2))*p1(1)+0.5*(1+gauss(2))*p2(1);
    ygp2 = 0.5*(1-gauss(2))*p1(2)+0.5*(1+gauss(2))*p2(2);
    xgp3 = 0.5*(1-gauss(3))*p1(1)+0.5*(1+gauss(3))*p2(1);
    ygp3 = 0.5*(1-gauss(3))*p1(2)+0.5*(1+gauss(3))*p2(2);
    
    % evaluate traction vector 'h' at gauss points
    hgp1 = [feval(traction_x,ygp1);feval(traction_y,xgp1)]; 
    hgp2 = [feval(traction_x,ygp2);feval(traction_y,xgp2)]; 
    hgp3 = [feval(traction_x,ygp3);feval(traction_y,xgp3)]; 
    
    % Due to linear triangular elements, every shape function along an
    % element edge is linear: N_1 = 0.5(1 - xsi), N_2 = 0.5(1 + xsi)
    % After evaluation at the gauss points 'gauss', they are the same for
    % every edge, if integrating along an element edge (Neumann boundary).
    % So, they are hard-coded here.
    N_11 = 0.5*(1-gauss(1));    % shape function 1 at gauss point 1
    N_12 = 0.5*(1+gauss(1));    % shape function 2 at gauss point 1
    N_21 = 0.5*(1-gauss(2));    % shape function 1 at gauss point 2
    N_22 = 0.5*(1+gauss(2));    % shape function 2 at gauss point 2
    N_31 = 0.5*(1-gauss(2));    % shape function 1 at gauss point 3
    N_32 = 0.5*(1+gauss(2));    % shape function 2 at gauss point 3
    
    
    % build shape function matrices 'Ngp1' and 'Ngp2' at gauss points 1 & 2
    Ngp1 = [N_11 0 N_12 0;
            0 N_11 0 N_12];
    Ngp2 = [N_21 0 N_22 0;
        0 N_21 0 N_22];
    Ngp3 = [N_31 0 N_32 0;
        0 N_31 0 N_32];
    
    force_values = Ngp1' * hgp1 * jac * weights(1) + ... 
        Ngp2' * hgp2 * jac * weights(2) + Ngp3' * hgp3 * jac * weights(3);
    force_id = [DOFs1(1) DOFs1(2) DOFs2(1) DOFs2(2)];
    
else
    error('MATLAB:XFEM:UnvalidForceShape', ...
        'No code for the traction distribution, you have chosen.');

end;