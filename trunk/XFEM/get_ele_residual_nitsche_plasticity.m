% get_ele_residual_nitsche_plasticity.m
%
% CALL: get_ele_residual_nitsche_plasticity(xcoords, ...
%         ycoords,seg_cut_info,endpoints,node, ...
%         x,y,dis,old_ndisp,id_dof,cutlist,maxngrains,totaldis, ...
%         id_eqns,GRAININFO_ARR,nodegrainmap,IFsliding_switch,alpha_t, ...
%         yieldstress,deltaload,dis_conv,old_ndisp_conv,IFsymmetrized)
%
% Evaluates the tangential residual contribution for an unsymmetric Nitsche
% formulation with perfect plasticity in tangential direction. The
% tangential traction is limited to a yield traction. This is done by a
% return-mapping algorithm.
%
% Input parameters:
%   xcoords             x-coordinates of current element
%   ycoords             y-coordinates of current element
%   seg_cut_info        information about current subsegment
%   endpoints           endpoints of current interface
%   node                connectivity between all nodes and elements
%   x                   x-coordinates of all nodes
%   y                   y-coordinates of all nodes
%   dis                 current re-assembled displacement vector
%   old_ndisp           copy of original current solution vector
%   id_dof              mapping between DOFs and enriching grains
%   cutlist             list of cut elements
%   maxngrains          number of all grains
%   fdisp               current solution vector
%   id_eqns             mappint between nodes and global DOFs (only for 
%                       this element)
%   GRAININFO_ARR       some data about the grains
%   nodegrainmap        mapping between nodes and grains
%   IFsliding_switch    indicates the sliding case
%   alpha _t            tangential stabilization parameter
%   yieldstress         yield stress / traction
%   deltaload           displacement increment for return mapping algorithm
%   dis_conv            total displacement vector of previous converged
%                       load step
%   old_ndisp_conv      unmodified total solution vector of previous 
%                       converged load step
%   IFsymmetrized       unsymmetric or symmetrized version of Nitsche's
%                       method with plasticity
%
% Returned variables:
%   res_nit_tangential  vector with tangential residual contribution
%   id                  id-array to prepare assembly
%   ttrac               tangential traction at the two gauss points
%   tgappl              tangential plastic gap at the two gauss points
%   ftrial              flow rule evaluated for trial state at the two
%                       gauss points

% Author: Matthias Mayr (08/2010)

function [res_nit_tangential id tangtrac tgappl f_trial] = ...
  get_ele_residual_nitsche_plasticity(xcoords, ...
    ycoords,seg_cut_info,endpoints,node, ...
    x,y,dis,old_ndisp,id_dof,cutlist,maxngrains,totaldis, ...
    id_eqns,GRAININFO_ARR,nodegrainmap,IFsliding_switch,alpha_t, ...
    yieldstress,deltaload,dis_conv,old_ndisp_conv,IFsymmetrized)

%% Initialize
elenodes = node(:,seg_cut_info.elemno);

res_nit_tang_base = zeros(6,1);       % residual for tangential constraints in base DOFs
res_nit_tang_enriched = zeros(12,1);  % residual for tangential constraints in enriched DOFs

enrich1 = zeros(1,3);             % Is there a first enrichment?
enrich2 = zeros(1,3);             % Is there a second enrichment?

intersection = seg_cut_info.xint; % intersection points of interface with 
                                  % element edges

normal = seg_cut_info.normal;     % vector normal to interface
tangent = seg_cut_info.tangent;   % vector tangential to interface

n_matrix = [normal(1) 0         normal(2);    % normal matrix to keep dimensions 
            0         normal(2) normal(1)];   %   consistent in discrete space
ntn = normal * normal';                       % "normal tensor normal" to project
                                              %   into normal direction

tgappl = [0 0];                   % plastic contribution to tangential gap
% plastictanggap = 0;

tangtrac = [0 0];                 % scalar values of current traction at gauss points
% ttrac_scalar = 0;

f_trial = [0 0];                  % flow rule of trial states at both gauss points
% f_trialgp = 0;
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
  enrich1(n) = id_dof(elenodes(n),3);

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
  enrich2(n) = id_dof(elenodes(n),5);  

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
%% GET DISCRETE CONSTITUTIVE MATRICES 'C1' and 'C2'
% grain 1
E = GRAININFO_ARR(pos_g).youngs;
pr = GRAININFO_ARR(pos_g).poisson;
fac = E/(1 - (pr)^2);
C1 = fac*[1.0,  pr,   0;
          pr,   1.0,  0.0;
          0,    0,    (1.0-pr)/2 ];

% grain 2
E = GRAININFO_ARR(neg_g).youngs;
pr = GRAININFO_ARR(neg_g).poisson;
fac = E/(1 - (pr)^2);
C2 = fac*[1.0,  pr,   0;
          pr,   1.0,  0.0;
          0,    0,    (1.0-pr)/2 ];
% ----------------------------------------------------------------------- %
%% COMPUTE AVERAGED STRESS
% Since the stress is constant in an element due to linear shape functions,
% this quantity can be computed outside the loop over the gauss points.

% Structure of 'stresse':
%   Dimension: 1x6xmaxngrains
%   Index 1:    'stresse' for one element (set to 1)
%   Index 2:    column 1    global element ID
%               column 2    x-coordinate of element centroid
%               column 3    y-coordinate of element centroid
%               column 4    xx-stress at element centroid
%               column 5    yy-stress at element centroid
%               column 6    xy-stress at element centroid
%   Index 3:    ID of grain, to which these values belong to (maxngrains =
%               maximum number of grains)

% Voigt-notation:
%   1. component: sigma_xx
%   2. component: sigma_yy
%   3. component: sigma_xy

% Due to the use of the Voigt notation, the product "sigma * normal" has
% to be computed with the modified matrix 'n_matrix' instead of the normal
% vector in order to preserve the right dimensionality.

% get stress in current element at the end of the previous converged load
% step
[~,stresse_conv] = post_process_better(node,x,y,seg_cut_info.elemno,dis_conv, ...
                old_ndisp_conv,id_dof,cutlist,maxngrains);
              
% get stress in current element with current displacement vector
[~,stresse_current] = post_process_better(node,x,y,seg_cut_info.elemno,dis, ...
                old_ndisp,id_dof,cutlist,maxngrains);
              
% compute stress increment
stresse = stresse_current - stresse_conv;
              
% stress increment in positive part of the element (Voigt-notation)
delta_stress1 = [ stresse(1,4,pos_g);
                  stresse(1,5,pos_g);
                  stresse(1,6,pos_g)];
          
          
% stress increment in negative part of the element (Voigt-notation)
delta_stress2 = [ stresse(1,4,neg_g);
            stresse(1,5,neg_g);
            stresse(1,6,neg_g)];

% compute averaged stress increment (Voigt-notation)
delta_stress_avg = 0.5 * (delta_stress1 + delta_stress2);%delta_stress1;%

% stress in positive part of the element (Voigt-notation)
stress1 = [ stresse_current(1,4,pos_g);
            stresse_current(1,5,pos_g);
            stresse_current(1,6,pos_g)];
          
          
% stress increment in negative part of the element (Voigt-notation)
stress2 = [ stresse_current(1,4,neg_g);
            stresse_current(1,5,neg_g);
            stresse_current(1,6,neg_g)];

% compute averaged stress (Voigt-notation)
stress_avg = 0.5 * (stress1 + stress2);%stress1;%
% ----------------------------------------------------------------------- %
%% GET AVERAGED 'CBhat' and 'CBtilde'
% Since the matrices 'Bhat' and 'Btilde' are constant in an element due to 
% linear shape functions, these quantities can be computed outside the loop
% over the gauss points.

% get B-matrices
[Bhat Btilde1 Btilde2] = getBmatrices(xcoords,ycoords, ...
  elenodes,nodegrainmap,pos_g,neg_g,pn_nodes,id_dof);


% Since 'Bhat1' = 'Bhat2', <CBhat> = <C>Bhat
CBhat_avg = 0.5 * (C1 + C2) * Bhat;%C1*Bhat;%

% Since 'Btilde1' ~= 'Btilde2', the average of 'Btilde' has to be
% considered, too.
CBtilde_avg = 0.5 * (C1 * Btilde1 + C2 * Btilde2);%C1*Btilde1;%
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
    N(1,2*b-1)    = N(1,2*b-1)    + Larea/Area;   % First enrichment
    N(2,2*b)      = N(2,2*b)      + Larea/Area;
    N(1,2*b-1+6)  = N(1,2*b-1+6)  + Larea/Area;   % Second enrichment
    N(2,2*b+6)    = N(2,2*b+6)    + Larea/Area;
  end;
  
  % multiply values belonging to nodes in the 'negative' grain with '-1'
  for c = 1:6
    N(1,2*c-1:2*c-1)  = N(1,2*c-1:2*c-1)*flg(c);
    N(2,2*c:2*c)      = N(2,2*c:2*c)    *flg(c);
  end;
  
  % evaluate total gap at current gauss point (elastic + plastic part)
  gap = evaluate_gap_gp(xn,yn,flg,xcoords,ycoords,totaldis', ...
    id_eqns,Area);
  
  % evaluate normal and tangential traction at gauss point
  % distinguish bewteen pure elastic and plastic cases
  if IFsliding_switch == 0 || IFsliding_switch == 1
    error('MATLAB:XFEM','This is the wrong routine.');
  elseif IFsliding_switch == 2
    % perfect plasticity
    [ttrac plastictanggap f_trialgp ttrac_scalar] = ...
      evaluate_traction_nitsche_gp_plastic(xn,yn,flg,xcoords,ycoords, ...
      id_eqns,Area,normal,tangent,alpha_t, ...
      seg_cut_info.tgapplconv(g),seg_cut_info.ttracconv(g),yieldstress, ...
      deltaload,delta_stress_avg,stress_avg,gap);
  else
    error('MATLAB:XFEM:UnvalidID','Unvalid sliding ID');
  end;
  
  % subtract the plastic part
  gap_el = gap - seg_cut_info.tgapplconv(g) * tangent;
  
  % tangential direction
  if IFsymmetrized == 0
    % unsymmetric
    res_nit_tang_base = res_nit_tang_base ...
      - CBhat_avg' * n_matrix' * (eye(2) - ntn) * gap_el * seg_jcob * weights(g);
    res_nit_tang_enriched = res_nit_tang_enriched ...
      - CBtilde_avg' * n_matrix' * (eye(2) - ntn) * gap_el * seg_jcob * weights(g);
  end;
  
  res_nit_tang_enriched = res_nit_tang_enriched + N' * ttrac * seg_jcob * weights(g);
  
  % store plastic contribution to tangential gap
  tgappl(g) = plastictanggap;
  
  % store plastic contribution to tangential gap
  tangtrac(g) = ttrac_scalar;
  
  % store flow rule of trial state
  f_trial(g) = f_trialgp;
end;

% assemble element residual from base and enriched contributions
res_nit_tangential = [res_nit_tang_base; 
                      res_nit_tang_enriched];
% ----------------------------------------------------------------------- %
%% Build id array
id = zeros(1,18);
id(1)  = id_eqns(1,1);  % original x dof
id(2)  = id_eqns(1,2);  % original y dof
id(3)  = id_eqns(2,1);  % original x dof
id(4)  = id_eqns(2,2);  % original y dof
id(5)  = id_eqns(3,1);  % original x dof
id(6)  = id_eqns(3,2);  % original y dof

id(7)  = id_eqns(1,3);  % 1st extra x dof
id(8)  = id_eqns(1,4);  % 1st extra y dof
id(9)  = id_eqns(2,3);  % 1st extra x dof
id(10) = id_eqns(2,4);  % 1st extra y dof
id(11) = id_eqns(3,3);  % 1st extra x dof
id(12) = id_eqns(3,4);  % 1st extra y dof

id(13) = id_eqns(1,5);  % 2nd extra x dof
id(14) = id_eqns(1,6);  % 2nd extra y dof
id(15) = id_eqns(2,5);  % 2nd extra x dof
id(16) = id_eqns(2,6);  % 2nd extra y dof
id(17) = id_eqns(3,5);  % 2nd extra x dof
id(18) = id_eqns(3,6);  % 2nd extra y dof
% ----------------------------------------------------------------------- %
end

