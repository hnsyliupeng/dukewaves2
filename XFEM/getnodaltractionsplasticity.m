% getnodaltractionsplasticity.m
%
% CALL: getnodaltractionsplasticity(xcoords,ycoords, ...
%   seg_cut_info,IFyieldstress,endpoints,id_dof,DOFs,NODEINFO_ARR)
%
% Compute nodal values of tangential tractions at the interface for 
% plasticity. An id-array for assembly will be provided, too.
%
% Input parameters:
%   xcoords           x-coordinates of the nodes
%   ycoords           y-coordinates of the nodes
%   seg_cut_info      informationt about current subsegment of interface
%   IFyieldstress     yield stress, given in input file
%   endpoints         points, that define the interface
%   id_dof            information, whether a node is enriched or not
%   DOFs              global DOFs, which belong to the three nodes
%   NODEINFO_ARR      some data about the nodes of the element
%
% Returned variables:
%   force_values      values of nodal forces for tangential traction
%   force_id          id-array with global DOFs for assembly into
%                     'big_force_traction'
%   tang_traction_max maximum tangential traction (limited due to
%                     plasticity)
%   he                physical length of current subsegment of interface
%

% Author: Matthias Mayr (06/2010)

function [force_values force_id tang_traction_max he] = ...
  getnodaltractionsplasticity(xcoords,ycoords,seg_cut_info, ...
  IFyieldstress,endpoints,id_dof,DOFs,NODEINFO_ARR)

% Initialize
xep = xcoords;
yep = ycoords;

enrich1 = zeros(1,3);
enrich2 = zeros(1,3);

% get intersection points
intersection = seg_cut_info.xint;

% get global ID of current element
eleID = seg_cut_info.elemno;

% Establish which nodes are "postively" enriched, and which reside in the 
% "negative" grain
pos_g = seg_cut_info.positive_grain;
neg_g = seg_cut_info.negative_grain;

[pn_nodes] = get_positive_new(eleID,pos_g,neg_g);

flg = [0 0 0 0 0 0];

% First enrichment
for n = 1:3     % loop over nodes
  % Get the "first" enrichment of node
  enrich1(n) = id_dof(n,3);

  if enrich1(n) == pos_g
    if (pn_nodes(n,1) == 1)
      flg(n) = 1;
    else
      flg(n) = 0;
    end;
  elseif enrich1(n) == neg_g
    if (pn_nodes(n,2) == 1)
      flg(n) = -1;
    else
      flg(n) = 0;
    end;        
  end;
end;

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
    end;
  elseif enrich2(n) == neg_g
    if (pn_nodes(n,2) == 1)
      flg(3 + n) = -1;
    else
      flg(3 + n) = 0;
    end;        
  end;
end;

% end points of intersection - direction doesn't matter - this is for the
% segment jacobian calculation
if all(size(intersection) == [2 2])     % no triple junction in the element
  % get the two intersection points
  p1 = intersection(1,:);
  p2 = intersection(2,:);
elseif all(size(intersection) == [1 2]) % triple junction in the element
  % get the two points, that define the subsegment
  p1 = intersection(1,:);

  % Second endpoint of segment is also end point of interface
  endpoint = endpoints(1,:);

  inside = polygon_contains_point_2d ( 3, [xep;yep], endpoint );

  if inside
    p2 = endpoint;
  else
    p2 = endpoints(2,:);
  end;
end;

% initialize
N = zeros(2,18);

% jacobian of segment to global
he = sqrt((p1(1)-p2(1))^2 + (p1(2)-p2(2))^2);
seg_jcob = he/2;

% Gauss points on segments
gauss = [-sqrt(3)/3 sqrt(3)/3];
weights = [1 1];

% The integration along the subsegment of the interface is done via Gauss'
% quadrature. Since the shape functions are linear and the tangential
% traction is constant in an element, you need two Gauss points for exact
% quadrature.

% loop over Gauss points to assemble N
for g = 1:2
  % Get real coordinates of gauss points
  xn = 0.5*(1-gauss(g))*p1(1)+0.5*(1+gauss(g))*p2(1);
  yn = 0.5*(1-gauss(g))*p1(2)+0.5*(1+gauss(g))*p2(2);

  for b = 1:3     % Evaluate shape functions
    % load node cooridates
    xes = xcoords;
    yes = ycoords;
    
    % Get coorindates of area opposite node of concern
    xes(b) = xn; 
    yes(b) = yn;

    Area = det([[1 1 1]' xep' yep'])/2;
    Larea = det([[1 1 1]' xes' yes'])/2;

    % Evaluate shape function for node 'b'
    N(1,2*b-1) = N(1,2*b-1) + Larea/Area*seg_jcob*weights(g);    % Base DOFs
    N(2,2*b)   = N(2,2*b)   + Larea/Area*seg_jcob*weights(g);
    N(1,2*b+5) = N(1,2*b+5) + Larea/Area*seg_jcob*weights(g);    % First enrichment
    N(2,2*b+6) = N(2,2*b+6) + Larea/Area*seg_jcob*weights(g);
    N(1,2*b+11) = N(1,2*b+11) + Larea/Area*seg_jcob*weights(g);  % Second enrichment
    N(2,2*b+12) = N(2,2*b+12) + Larea/Area*seg_jcob*weights(g);
  end;
end;

% set some values to zero depending on, whether they are "positively" or
% "negatively" enriched.
for c = 1:6
  N(:,6 + 2*c-1:6+2*c) = N(:,6 + 2*c-1:6+2*c)*flg(c);
end;
% here N contains the evatluated shape functions for two possible
% enrichments. If there is only one enrichment in this element, the
% corresponding entries in N are set to zero by the latter for-loop.

% set the base DOFs to zero, since the traction depends only on the
% enriched DOFs
for c=1:6
  N(1:2,c) = 0;
end;

% vector of tangential traction (computed via yield stress)
tang_traction_max = IFyieldstress * seg_cut_info.tangent * he;

% compute nodal force values
force_values = (-1) *  N' * tang_traction_max;

% Set all nodal forces to zero, whose nodes don't reside in the enriching 
% grain
for i=1:3
  if id_dof(i,3) ~= NODEINFO_ARR(1,i).grain
    force_values([2*i 2*i-1 2*i+5 2*i+6 2*i+11 2*i+12]) = 0;
  end;
end;

% set up the id-array
force_id = [DOFs(1,1:2) DOFs(2,1:2) DOFs(3,1:2) ... % base DOFs
  DOFs(1,3:4) DOFs(2,3:4) DOFs(3,3:4) ...           % first enrichment
  DOFs(1,5:6) DOFs(2,5:6) DOFs(3,5:6)];             % second enrichment  

end