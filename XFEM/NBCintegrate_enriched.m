% NBCintegrate_enriched.m
%
% CALL: NBCintegrate_enriched(BOUNDARY,FORCE,DOFs1,DOFs2,seg_cut_info,NODEINFO_ARR,
%       id_dof)
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
%   seg_cut_info      structure with information about interface and
%                     enrichment
%   NODEINFO_ARR      structure with information about the two nodes
%   id_dof            global degrees of freedom of enriched nodes
%
% Returned variables
%   force_values      element load vector (nodal forces)
%   force_id          ID array for assembly

% Author: Matthias Mayr (05/2010)
function [force_values,force_id] = NBCintegrate_enriched(BOUNDARY,FORCE,...
  DOFs1,DOFs2,seg_cut_info,NODEINFO_ARR,id_dof)

if strcmpi(FORCE.shape,'constant')
  % Now, use 2-point-gauss-integration to evaluate the integral over the
  % Neumman boundary: \int _{\Gamma _h} w _i h _i d \Gamma
  % with  w   test function
  %       h   traction vector
  
  % define some values for 2-point-gauss integration
  gauss = [-sqrt(3)/3 sqrt(3)/3];
  weights = [1.0 1.0];

  % compute length of Neumann boundary
  len = norm(FORCE.coords(1,:) - FORCE.coords(2,:));
  
  % get coords of element nodes
  p1 = BOUNDARY.coords(1,:);
  p2 = BOUNDARY.coords(2,:);
  
  % get coords of intersection point
  p = seg_cut_info.xint;
  for i=1:2
    dist(i) = line_exp_point_dist_2d(p1,p2,p);
  end;
  if dist(1)<dist(2)  % choose point, which is closer to element edge
    p = p(1,:);
  else
    p = p(2,:);
  end;
  
  % jacobian depends on, which part of the element edge is in the enriched
  % grain.
  % compute jacobian
  if NODEINFO_ARR(1).grain == id_dof(1,3)
    jac = norm(p1 - p)/2;
  else
    jac = norm(p2 - p)/2;
  end;
    
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
  Ngp1 = [N_11 0 N_12 0;    % elementary matrix
          0 N_11 0 N_12];
  Ngp2 = [N_21 0 N_22 0;    % elementary matrix
      0 N_21 0 N_22];
%   % matrices for base and enriched DOFs
%   Ngp1 = [Ngp1(:,1:2) Ngp1(:,1:2) Ngp1(:,3:4) Ngp1(:,3:4)];       
%   Ngp2 = [Ngp2(:,1:2) Ngp2(:,1:2) Ngp2(:,3:4) Ngp2(:,3:4)];
  
  % set some values in 'Ngp1' and 'Ngp2' to zero, if the corresponding node
  % is not enriched
  % evaluated at gauss point 1
  if DOFs1(3) == 0, Ngp1(1,1)=0; end; % node 1 --> x
  if DOFs1(4) == 0, Ngp1(2,2)=0; end; % node 1 --> y
  if DOFs2(3) == 0, Ngp1(1,3)=0; end; % node 2 --> x
  if DOFs2(4) == 0, Ngp1(2,4)=0; end; % node 2 --> y
  % evaluated at gauss point 2
  if DOFs1(3) == 0, Ngp2(1,1)=0; end; % node 1 --> x
  if DOFs1(4) == 0, Ngp2(2,2)=0; end; % node 1 --> y
  if DOFs2(3) == 0, Ngp2(1,3)=0; end; % node 2 --> x
  if DOFs2(4) == 0, Ngp2(2,4)=0; end; % node 2 --> y
  
  % traction vector, evaluated at both gauss points
  hgp1 = [FORCE.values(1); FORCE.values(2)] ./ len;
  hgp2 = [FORCE.values(1); FORCE.values(2)] ./ len;
  
  % compute element load vector and index array for assembly
  force_values = Ngp1' * hgp1 * jac * weights(1) + ... 
    Ngp2' * hgp2 * jac * weights(2);
  force_id = [DOFs1(3) DOFs1(4) DOFs2(3) DOFs2(4)];
  
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
  traction_x = @(x) mx*x + tx;
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

  % compute element load vector and index array for assembly
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
  jac = norm(BOUNDARY.coords(1,:) - BOUNDARY.coords(2,:)) / 2;


  % comput length of Neumann boundary
  len = norm(FORCE.coords(1,:) - FORCE.coords(2,:));
  lenhalf = len / 2; % equals to 'h' from description above

  % construct function, that describes Lagrange polynomial for parabolic
  % force distribution
  fx = FORCE.values(1,1);    % value of x-force at vertex
  fy = FORCE.values(1,2);    % value of y-force at vertex

  % define functions for prabolic traction distribution
  % x-traction as a function of x (set to zero)
  traction_x = @(x) 0;

  % y-traction as a function of y
  traction_y = @(y) (-1) * fy / (lenhalf^2) * (y^2) + fy;

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
  hgp1 = [feval(traction_x,xgp1);feval(traction_y,ygp1)]; 
  hgp2 = [feval(traction_x,xgp2);feval(traction_y,ygp2)]; 
  hgp3 = [feval(traction_x,xgp3);feval(traction_y,ygp3)]; 

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