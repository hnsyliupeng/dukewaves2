% get_lag_mults_for_penalty_alternative.m
%
% CALL: get_lag_mults_for_penalty_alternative(xcoords, ...
%         ycoords,seg_cut_info,endpoints,id_dof,id_eqns,fdisp)
%
% This method computes the traction at the interface via the penalty method
% using the one-integral formulation.
%
% Input arguments:
%   xcoords             x-coordinates of element's nodes
%   ycoords             y-coordinates of element's nodes
%   seg_cut_info        information about the interface in the element
%   endpoints           endpoints of interface
%   id_dof              shows, if a node is enriched or not
%   id_eqns             mapping between element's nodes and global DOFs
%   fdisp               displacements in base and extra DOFs of the element
%
% Returned variables
%   lagmult         Vector of Lagrange multipliers
%

% Author: Matthias Mayr (04/2010)

function [normaltraction tangentialtraction] = ...
  get_lag_mults_for_penalty_alternative(xcoords, ...
  ycoords,seg_cut_info,endpoints,id_dof,id_eqns,fdisp)

% Initialize
enrich1 = zeros(1,3);
enrich2 = zeros(1,3);

N = zeros(2,12);  % shape function matrix

intersection = seg_cut_info.xint;

% Establish a set of flags
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

  inside = polygon_contains_point_2d ( 3, [xcoord;ycoord], endpoint );

  if inside
    p2 = endpoint;
  else
    p2 = endpoints(2,:);
  end
end

% Get real coordinates midpoint of subsegment
xn = 0.5 * (p1(1) + p2(1));
yn = 0.5 * (p1(2) + p2(2));

% get area of element
Area = det([[1 1 1]' xcoords' ycoords'])/2;

% Evaluate shape functions
for b = 1:3     
  % get area opposite of node of concern
  xes = xcoords;
  yes = ycoords;
  xes(b) = xn; 
  yes(b) = yn;
  Larea = det([[1 1 1]' xes' yes'])/2;

  % Evaluate shape function for node 'b'
  N(1,2*b-1) = N(1,2*b-1) + Larea/Area;    % First enrichment
  N(2,2*b)   = N(2,2*b)   + Larea/Area;
  N(1,2*b+5) = N(1,2*b+5) + Larea/Area;    % Second enrichment
  N(2,2*b+6) = N(2,2*b+6) + Larea/Area;
end

% set some values to zero depending on, whether they are "positively" or
% "negatively" enriched.
for c = 1:6
  N(:,2*c-1:2*c) = N(:,2*c-1:2*c)*flg(c);
end
% here N contains the evatluated shape functions for two possible
% enrichments. If there is only one enrichment in this element, the
% corresponding entries in N are set to zero by the latter for-loop.

%--------------------------------------------------------------------------

% To decide, if you need one or two enrichments, you cannot look at the
% element, you have to look at each node of the element, because the node
% can be enriched two times, also if this element is only cut by one
% interface. This is the fact, if the double-enriched node belongs to a
% neighboured element, that contains a triple junction.
%
% Check every node, if it is enriched or not. If it is enriched only once, 
% add two entries to 'localdis' and fill up with two zeros. If it is 
% enriched twice, fill in the extra DOFs.

% Check first node
index1 = id_eqns(1,3:6);
if all(index1)      % node 1 of element 'parent' is enriched twice
  localdis1 = [fdisp(index1(1:2)); zeros(4,1); fdisp(index1(3:4)); zeros(4,1)];
else                % node 1 of element 'parent' is enriched once
  localdis1 = [fdisp(index1(1:2)); zeros(10,1)];
end;

% Check second node
index2 = id_eqns(2,3:6);
if all(index2)      % node 2 of element 'parent' is enriched twice
  localdis2 = [zeros(2,1); fdisp(index2(1:2)); zeros(4,1); ...
    fdisp(index2(3:4)); zeros(2,1)];
else                % node 2 of element 'parent' is enriched once
  localdis2 = [zeros(2,1); fdisp(index2(1:2)); zeros(8,1)];
end;

% Check third node
index3 = id_eqns(3,3:6);
if all(index3)      % node 3 of element 'parent' is enriched twice
  localdis3 = [zeros(4,1); fdisp(index3(1:2)); zeros(4,1); ...
    fdisp(index3(3:4))];
else                % node 3 of element 'parent' is enriched once
  localdis3 = [zeros(4,1); fdisp(index3(1:2)); zeros(6,1)];
end;

% Assemble 'localdis'
localdis = localdis1' + localdis2' + localdis3';

% compute traction
traction = (N * localdis')';

normaltraction = traction * seg_cut_info.normal;
tangentialtraction = traction * seg_cut_info.tangent;

end

