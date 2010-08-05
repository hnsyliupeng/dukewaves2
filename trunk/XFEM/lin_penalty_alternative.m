%% lin_penalty_alternative.m
%
% CALL: lin_penalty_alternative(xcoords,ycoords,id_eqns,id_dof, ...
%         seg_cut_info,endpoints)
%
% Computes the penalty-contributions for constraints in normal and
% tangential direction to the tangent matrix 'tangentmatrix'. In opposite
% to 'lin_penalty.m', a single integral formulation is used which requires
% the following scaling of the parameter 'alpha'
%   1/h^2               for penalty method
%   1/h                 for Nitsche's method
% An alternative form with two integrals is implemented in 'lin_penalty.m'.
%
% Input arguments:
%   xcoords             x-coordinates of element's nodes
%   ycoords             y-coordinates of element's nodes
%   parent              global element ID of current element
%   id_eqns             mapping between element's nodes and global DOFs
%   id_dof              shows, if a node is enriched or not
%   seg_cut_info        information about the interface in the element
%   endpoints           endpoints of interface
%   alpha_n             normal penalty parameter
%   alpha_t             tangential penalty parameter
%
% Returned variables
%   pen_normal          matrix for penalty constraints in normal direction
%   pen_tangent         matrix for penalty constraints in tangential direction
%   id                  id-array to enable assembly into 'tangentmatrix'
%

% Author: Matthias Mayr (07/2010)

function [pen_normal pen_tangent id] =... 
    lin_penalty_alternative(xcoords,ycoords,id_eqns,id_dof, ...
    seg_cut_info,endpoints,alpha_n,alpha_t)
%% Initialize
% penalty stiffness contribution
pen_normal = zeros(12,12);             % matrix for normal constraints
pen_tangent = zeros(12,12);            % matrix for tangential constraints

enrich1 = zeros(1,3);               % Is there a first enrichment?
enrich2 = zeros(1,3);               % Is there a second enrichment?

normal = seg_cut_info.normal;       % vector normal to interface

intersection = seg_cut_info.xint;   % intersection points between interface
                                    % and element edges
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

% 2 Gauss points on segments
gauss = [-sqrt(3)/3 sqrt(3)/3];
weights = [1 1];

% get area of element
Area = det([[1 1 1]' xcoords' ycoords'])/2;
% ----------------------------------------------------------------------- %
%% loop over Gauss points to assemble N
for g = 1:length(gauss)
  % Get real coordinates of gauss points
  xn = 0.5*(1-gauss(g))*p1(1)+0.5*(1+gauss(g))*p2(1);
  yn = 0.5*(1-gauss(g))*p1(2)+0.5*(1+gauss(g))*p2(2);

  % reset shape function matrix 
  N = zeros(2,12);
  
  % Evaluate shape functions
  for b = 1:3     
    % Get coorindates of area opposite node of concern
    xes = xcoords;
    yes = ycoords;
    xes(b) = xn; 
    yes(b) = yn;
    Larea = det([[1 1 1]' xes' yes'])/2;

    % Evaluate shape function
    N(1,2*b-1) = N(1,2*b-1) + Larea/Area;    % First enrichment
    N(2,2*b)   = N(2,2*b)   + Larea/Area;
    N(1,2*b+5) = N(1,2*b+5) + Larea/Area;    % Second enrichment
    N(2,2*b+6) = N(2,2*b+6) + Larea/Area;
  end
  
  % multiply values belonging to nodes in the 'negative' grain with '-1'
  for c = 1:6
      N(:,2*c-1:2*c) = N(:,2*c-1:2*c)*flg(c);
  end
  
  % assemble normal penalty matrix
  pen_normal = pen_normal ...
    + alpha_n * N' * (normal * normal') * N * weights(g);
  
  % assemble tangential penalty matrix
  if seg_cut_info.f_trial(g) <= 0   % only, if there is a tangential 
                                  % stiffness at current gauss point
    pen_tangent = pen_tangent ...
      + alpha_t * N' * (eye(2) - normal * normal') * N * weights(g);
  end;
end;

% multiply with jacobian due to gauss quadrature
pen_normal = pen_normal .* seg_jcob;
pen_tangent = pen_tangent .* seg_jcob;
% ----------------------------------------------------------------------- %      
%% Build id array
% get DOFs for first enrichment
id(1) = id_eqns(1,3);  % 1st extra x dof
id(2) = id_eqns(1,4);  % 1st extra y dof
id(3) = id_eqns(2,3);  % 1st extra x dof
id(4) = id_eqns(2,4);  % 1st extra y dof
id(5) = id_eqns(3,3);  % 1st extra x dof
id(6) = id_eqns(3,4);  % 1st extra y dof

% get DOFs for second enrichment
id(7)  = id_eqns(1,5);  % 2nd extra x dof
id(8)  = id_eqns(1,6);  % 2nd extra y dof
id(9)  = id_eqns(2,5);  % 2nd extra x dof
id(10) = id_eqns(2,6);  % 2nd extra y dof
id(11) = id_eqns(3,5);  % 2nd extra x dof
id(12) = id_eqns(3,6);  % 2nd extra y dof
