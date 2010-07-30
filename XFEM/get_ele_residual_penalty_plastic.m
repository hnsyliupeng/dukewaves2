% get_ele_residual_penalty_plastic.m
%
% CALL: get_ele_residual_penalty_plastic()
%
% Computes the nodal values for constant contribution to the residual for a
% penalty method. It is neded as the plastic traction if a elasto-plastic
% interfacial law is used.
%
% Input parameters:
%   xcoords             x-coordinates of element's nodes
%   ycoords             y-coordinates of element's nodes
%   seg_cut_info        information about the interface in the element
%   endpoints           endpoints of interface
%   id_dof              shows, if a node is enriched or not
%   id_eqns             mapping between element's nodes and global DOFs
%   totaldis            displacements in base and extra DOFs of the element
%
% Returned variables:
%
%

% Author: Matthias Mayr (07/2010)

function [constant_tangent] = get_ele_residual_penalty_plastic(xcoords, ...
  ycoords,seg_cut_info,endpoints,id_dof,id_eqns,fdisp)
%% Initialize

intersection = seg_cut_info.xint; % intersection points of interface with 
                                  % element edges

normal = seg_cut_info.normal;     % vector normal to interface

% shape function matrix for first and second enrichments
N = zeros(2,12);
% ----------------------------------------------------------------------- %
%% Establish a set of flags
flg = [0 0 0 0 0 0];

% Establish which nodes are "postively" enriched, and
% which reside in the "negative" grain
pos_g = seg_cut_info.positive_grain;
neg_g = seg_cut_info.negative_grain;

[pn_nodes] =... 
   get_positive_new(seg_cut_info.elemno,pos_g,neg_g);

% First enrichment
for n = 1:3     % loop over nodes
  % Get the "first" enrichment of node
  enrich1(n) = id_dof(n,3);

  if enrich1(n) == pos_g
    if (pn_nodes(n,1) == 1)
      flg(n) = 1;
    else
      flg(n) = 0;
    end
  elseif enrich1(n) == neg_g
    if (pn_nodes(n,2) == 1)
      flg(n) = -1;
    else
      flg(n) = 0;
    end        
  end
end

% Second Enrichment
for n = 1:3     % loop over nodes
  % Get the "second" enrichment of nodes
  enrich2(n) = id_dof(n,5);  

  if enrich2(n) == pos_g  % If this enrichment corresponds 
                          % to the positive grain
    if (pn_nodes(n,1) == 1)
      flg(3 + n) = 1;
    else
      flg(3 + n) = 0;
    end

  elseif enrich2(n) == neg_g
    if (pn_nodes(n,2) == 1)
      flg(3 + n) = -1;
    else
      flg(3 + n) = 0;
    end        
  end
end
% ----------------------------------------------------------------------- %
%% PREPARE GAUSS QUADRATURE
% end points of intersection - direction doesn't matter - this is for the
% segment jacobian calculation

if all(size(intersection) == [2 2])
  p1 = intersection(1,:);
  p2 = intersection(2,:);
elseif all(size(intersection) == [1 2])
  p1 = intersection(1,:);

  % Second endpoint of segment is also end point of interface
  endpoint = endpoints(1,:);

  inside = polygon_contains_point_2d ( 3, [xcoords;ycoords], endpoint );

  if inside
    p2 = endpoint;
  else
    p2 = endpoints(2,:);
  end
end

% jacobian of segment to global
he = sqrt((p1(1)-p2(1))^2 + (p1(2)-p2(2))^2);
seg_jcob = he/2;

% 2 point Gauss quadrature
gauss = [-sqrt(3)/3 sqrt(3)/3];
weights = [1 1];

% get area of element
Area = det([[1 1 1]' xcoords' ycoords'])/2;
% ----------------------------------------------------------------------- %
%% loop over Gauss points to assemble 'N'
for g = 1:length(gauss)
  %% assemble N
  % Get real coordinates of gauss points
  xn = 0.5*(1-gauss(g))*p1(1)+0.5*(1+gauss(g))*p2(1);
  yn = 0.5*(1-gauss(g))*p1(2)+0.5*(1+gauss(g))*p2(2);

  for b = 1:3     % Evaluate shape functions
    % Get coorindates of area opposite node of concern
    xes = xcoords;
    yes = ycoords;
    xes(b) = xn; 
    yes(b) = yn;
    Larea = det([[1 1 1]' xes' yes'])/2;

    % Evaluate shape function
    N(1,2*b-1) = N(1,2*b-1) + Larea/Area .* weights(g) * seg_jcob;   % First enrichment
    N(2,2*b)   = N(2,2*b)   + Larea/Area .* weights(g) * seg_jcob;
    N(1,2*b+5) = N(1,2*b+5) + Larea/Area .* weights(g) * seg_jcob;   % Second enrichment
    N(2,2*b+6) = N(2,2*b+6) + Larea/Area .* weights(g) * seg_jcob;
  end;
end;
  
% multiply values belonging to nodes in the 'negative' grain with '-1'
for c = 1:6
  N(:,2*c-1:2*c) = N(:,2*c-1:2*c)*flg(c);
end;
  
% Check first node
index1 = id_eqns(1,3:6);
if all(index1)      % node 1 of element 'parent' is enriched twice
  localdis1 = [fdisp(index1(1:2)) zeros(1,4) fdisp(index1(3:4)) zeros(1,4)];
else                % node 1 of element 'parent' is enriched once
  localdis1 = [fdisp(index1(1:2)) zeros(1,10)];
end;

% Check second node
index2 = id_eqns(2,3:6);
if all(index2)      % node 2 of element 'parent' is enriched twice
  localdis2 = [zeros(1,2) fdisp(index2(1:2)) zeros(1,4) ...
    fdisp(index2(3:4)) zeros(1,2)];
else                % node 2 of element 'parent' is enriched once
  localdis2 = [zeros(1,2) fdisp(index2(1:2)) zeros(1,8)];
end;

% Check third node
index3 = id_eqns(3,3:6);
if all(index3)      % node 3 of element 'parent' is enriched twice
  localdis3 = [zeros(1,4) fdisp(index3(1:2)) zeros(1,4) ...
    fdisp(index3(3:4))];
else                % node 3 of element 'parent' is enriched once
  localdis3 = [zeros(1,4) fdisp(index3(1:2)) zeros(1,6)];
end;

% Assemble 'localdis'
localdis = localdis1' + localdis2' + localdis3';

% compute tangent vector
constant_tangent = N' * (eye(2) - normal * normal') * N * localdis;
end

