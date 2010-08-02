% evaluate_traction_penalty_gp_plastic.m
%
% CALL: evaluate_gap_penalty_gp_plastic(xn,yn,flg, ...
%         xcoords,ycoords,fdisp,id_eqns,Area,normal,tangent,alpha_n, ...
%         alpha_t,tgapplconv,ttracconv,yieldstress,deltaload)
%
% This method computes the Lagrange multipliers (internal forces at the
% interface) via the penalty method. They are obtained by integrating over
% the gap, so they will be piecewise constant.
%
% Input arguments:
%   xn              x-coordinate of current gauss point
%   yn              y-coordinate of current gauss point
%   xcoords         x-coordinates of current element's nodes
%   ycoords         y-coordinates of current element's nodes
%   fdisp           nodal displacements of current element
%   id_eqns         mapping between nodes and DOFs for current element
%   Area            area of the current element
%   normal          vector normal to the interface
%   tangent         vector tangential to the interface
%   alpha_n         normal penalty parameter
%   alpha_t         tangential penalty parameter
%   tgapplconv      converged tangential plastic gap of previous load step
%   ttracconv       converged tangential traction of previous load step
%   yieldstress     yield stress
%
% Returned variables
%   ntrac           normal traction
%   ttrac           tangential traction
%   tgapppl         plastic contribution to tangential gap
%   f_trial         flow rule evaluated at trial state (assumed elasticity)
%   ttrac_scalar    scalar value of tangential tractio
%

% Author: Matthias Mayr (07/2010)

function [ntrac ttrac tgappl f_trial ttrac_scalar] = ...
  evaluate_traction_penalty_gp_plastic(xn,yn,flg,xcoords,ycoords,fdisp, ...
  id_eqns,Area,normal,tangent,alpha_n,alpha_t,tgapplconv,ttracconv, ...
  yieldstress,deltaload)

N = zeros(2,12);

for b = 1:3     % Evaluate shape functions
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

%% evaluate gap
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

% compute gap
gap = N * localdis';
% ----------------------------------------------------------------------- %
%% evaluate 'tgapcurrent'
% Check first node
index1 = id_eqns(1,3:6);
if all(index1)      % node 1 of element 'parent' is enriched twice
  localdis1 = [deltaload(index1(1:2)); zeros(4,1); deltaload(index1(3:4)); zeros(4,1)];
else                % node 1 of element 'parent' is enriched once
  localdis1 = [deltaload(index1(1:2)); zeros(10,1)];
end;

% Check second node
index2 = id_eqns(2,3:6);
if all(index2)      % node 2 of element 'parent' is enriched twice
  localdis2 = [zeros(2,1); deltaload(index2(1:2)); zeros(4,1); ...
    deltaload(index2(3:4)); zeros(2,1)];
else                % node 2 of element 'parent' is enriched once
  localdis2 = [zeros(2,1); deltaload(index2(1:2)); zeros(8,1)];
end;

% Check third node
index3 = id_eqns(3,3:6);
if all(index3)      % node 3 of element 'parent' is enriched twice
  localdis3 = [zeros(4,1); deltaload(index3(1:2)); zeros(4,1); ...
    deltaload(index3(3:4))];
else                % node 3 of element 'parent' is enriched once
  localdis3 = [zeros(4,1); deltaload(index3(1:2)); zeros(6,1)];
end;

% Assemble 'localdis'
localdis = localdis1' + localdis2' + localdis3';

% compute 'tgapcurrent'
tgapcurrent = tangent' * N * localdis';
% ----------------------------------------------------------------------- %
%% compute tractions
% compute normal traction
ntrac = alpha_n * (normal * normal') * gap;

% compute tangential traction, using a return map algorithm
[ttrac_scalar tgappl f_trial] = returnmappingpenalty(alpha_t,tgapplconv, ...
  ttracconv,yieldstress,tgapcurrent);
ttrac = ttrac_scalar * tangent;
% ttrac_vec = (eye(2) - normal * normal') * gap;
% ttrac = ttrac_scalar / norm(ttrac_vec) * ttrac_vec;
% ----------------------------------------------------------------------- %
end

