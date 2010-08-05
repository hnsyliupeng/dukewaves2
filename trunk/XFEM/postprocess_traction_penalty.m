% postprocess_traction_penalty.m
%
% CALL: postprocess_traction_penalty()
%
% Evaluates the normal tractions at the two gauss points in an subsegment
% of the interface for a penalty method.
%
% Input arguments:
%   xcoords             x-coordinates of element's nodes
%   ycoords             y-coordinates of element's nodes
%   seg_cut_info        information about the interface in the element
%   endpoints           endpoints of interface
%   id_dof              shows, if a node is enriched or not
%   id_eqns             mapping between element's nodes and global DOFs
%   fdisp               displacements in base and extra DOFs of the element
%   alpha_n             normal penalty parameter
%   alpha_t             tangential penalty parameter
%
% Returned variables
%   ntrac
%

% Author: Matthias Mayr (08/2010)

function [ntrac ttrac] = ...
  postprocess_traction_penalty(xcoords,ycoords,seg_cut_info, ...
  id_dof,id_eqns,fdisp,alpha_n,alpha_t,endpoints)

%% INITIALIZE
enrich1 = zeros(1,3);
enrich2 = zeros(1,3);

intersection = seg_cut_info.xint;

normal = seg_cut_info.normal;
tangent = seg_cut_info.tangent;

ntrac = [0 0];
ttrac = [0 0];

%% ESTABLISH A SET OF FLAGS
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
%% get nodal displacement vector
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
% ----------------------------------------------------------------------- %
%% PREPARE FINDING THE GAUSS POINTS
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

  inside = polygon_contains_point_2d ( 3, [xcoords;ycoords], endpoint );

  if inside
    p2 = endpoint;
  else
    p2 = endpoints(2,:);
  end
end

% Gauss points on segments
gauss = [-sqrt(3)/3 sqrt(3)/3];

% compute area of element
Area = det([[1 1 1]' xcoords' ycoords'])/2;

%% LOOP OVER GAUSS POINTS
for g = 1:2
  % reset shape function matrix
  N = zeros(2,12);
  
  % Get real coordinates of gauss point 'g'
  xn = 0.5*(1-gauss(g))*p1(1)+0.5*(1+gauss(g))*p2(1);
  yn = 0.5*(1-gauss(g))*p1(2)+0.5*(1+gauss(g))*p2(2);
  
  % Evaluate shape functions and assemble 'N'
  for b = 1:3     
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
   
  % compute scalar values of tractions at gauss point 'g'
  ntrac(g) = alpha_n .* normal' * N * localdis';
  ttrac(g) = alpha_t .* tangent' * N * localdis';
end

end

