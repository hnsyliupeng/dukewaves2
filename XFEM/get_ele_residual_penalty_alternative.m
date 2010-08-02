% get_ele_residual_penalty_alternative.m
%
% CALL: get_ele_residual_penalty_alternative(xcoords,ycoords, ...
%         seg_cut_info,endpoints,id_dof,id_eqns,totaldis)
%
% Computes the contribution of the contraints to the residual, depending on
% the method of constraint enforcement.
% 
% Input arguments:
%   xcoords             x-coordinates of element's nodes
%   ycoords             y-coordinates of element's nodes
%   seg_cut_info        information about the interface in the element
%   endpoints           endpoints of interface
%   id_dof              shows, if a node is enriched or not
%   id_eqns             mapping between element's nodes and global DOFs
%   totaldis            displacements in base and extra DOFs of the element
%   alpha_n             normal penalty parameter
%   alpha_t             tangential penalty parameter
%   yieldstress         yield stress
%   deltaload           
%   sliding_switch      indicates the sliding mode
%
% Returned variables
%   penalty_normal      penalty residual in normal direction
%   penalty_tangent     penalty residual in tangential direction 
%   id                  id-array to enable assembly into 'residual'
%   tgappl              plastic contribution to tangential gap
%   tangtrac            tangential tractions at each gauss point
%   f_trial             evalutated flow rule of trial state at both gauss 
%                       points
%

% Author: Matthias Mayr (07/2010)

function [penalty_normal penalty_tangent id tgappl tangtrac f_trial] = ...
  get_ele_residual_penalty_alternative(xcoords,ycoords,seg_cut_info, ...
  endpoints,id_dof,id_eqns,totaldis,alpha_n,alpha_t,yieldstress, ...
  deltaload,sliding_switch)
%% Initialize
penalty_normal = zeros(12,1);    % matrix for normal constraints
penalty_tangent = zeros(12,1);   % matrix for tangential constraints

enrich1 = zeros(1,3);             % Is there a first enrichment?
enrich2 = zeros(1,3);             % Is there a second enrichment?

intersection = seg_cut_info.xint; % intersection points of interface with 
                                  % element edges

normal = seg_cut_info.normal;     % vector normal to interface
tangent = seg_cut_info.tangent;   % vector tangential to interface

tgappl = [0 0];                   % plastic contribution to tangential gap
plastictanggap = 0;

tangtrac = [0 0];                 % scalar values of current traction at gauss points
ttrac_scalar = 0;

f_trial = [0 0];                  % flow rule of trial states at both gauss points
f_trialgp = 0;
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
%% PREPARE GAUSS QUADRATURE (TWO GAUSS POINTS)
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

% 2 gauss points on segments, since the product of linear shape function
% matrices is integrated, i.e. the integrand is quadratic, so at least 2
% gauss points are required for an exact integration.
gauss = [-sqrt(3)/3 sqrt(3)/3];
weights = [1 1];

% get area of element
Area = det([[1 1 1]' xcoords' ycoords'])/2;
% ----------------------------------------------------------------------- %
%% loop over Gauss points to assemble 'N'
for g = 1:length(gauss)
  % Get real coordinates of gauss points
  xn = 0.5*(1-gauss(g))*p1(1)+0.5*(1+gauss(g))*p2(1);
  yn = 0.5*(1-gauss(g))*p1(2)+0.5*(1+gauss(g))*p2(2);

  % reset shape function matrix 'N'
  N = zeros(2,12);
  
  % Evaluate shape functions and assemble 'N'
  for b = 1:3     
    % Get coorindates of area opposite node of concern
    xes = xcoords;
    yes = ycoords;
    xes(b) = xn; 
    yes(b) = yn;
    Larea = det([[1 1 1]' xes' yes'])/2;

    % Evaluate shape function
    N(1,2*b-1) = N(1,2*b-1) + Larea/Area;   % First enrichment
    N(2,2*b)   = N(2,2*b)   + Larea/Area;
    N(1,2*b+5) = N(1,2*b+5) + Larea/Area;   % Second enrichment
    N(2,2*b+6) = N(2,2*b+6) + Larea/Area;
  end;
  
  % multiply values belonging to nodes in the 'negative' grain with '-1'
  for c = 1:6
    N(:,2*c-1:2*c) = N(:,2*c-1:2*c)*flg(c);
  end;
  
  % evaluate normal and tangential traction at gauss point
  % distinguish bewteen pure elastic and plastic cases
  if sliding_switch == 0 || sliding_switch == 1
    [ntrac ttrac] = evaluate_traction_penalty_gp(xn,yn,flg,xcoords, ...
      ycoords,totaldis,id_eqns,Area,normal,alpha_n,alpha_t);
  elseif sliding_switch == 2
    % perfect plasticity
    [ntrac ttrac plastictanggap f_trialgp ttrac_scalar] = ...
      evaluate_traction_penalty_gp_plastic(xn,yn,flg,xcoords,ycoords, ...
      totaldis,id_eqns,Area,normal,tangent,alpha_n,alpha_t, ...
      seg_cut_info.tgapplconv(g),seg_cut_info.ttracconv(g),yieldstress,deltaload);
  else
    error('MATLAB:XFEM:UnvalidID','Unvalid sliding ID');
  end;
  % normal direction
  penalty_normal  = penalty_normal + N' * ntrac * weights(g);

  % tangential direction
  penalty_tangent = penalty_tangent + N' * ttrac * weights(g);
  
  % store plastic contribution to tangential gap
  tgappl(g) = plastictanggap;
  
  % store plastic contribution to tangential gap
  tangtrac(g) = ttrac_scalar;
  
  % store flow rule of trial state
  f_trial(g) = f_trialgp;
end;

% multiply with jacobian due to gauss quadrature
penalty_normal = penalty_normal .* seg_jcob;
penalty_tangent = penalty_tangent .* seg_jcob;
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
% ----------------------------------------------------------------------- %
end

