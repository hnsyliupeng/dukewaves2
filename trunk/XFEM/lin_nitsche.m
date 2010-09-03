% lin_nitsche.m
%
% CALL: lin_nitsche(xcoords,ycoords,seg_cut_info,endpoints,node,id_dof, ...
%         id_eqns,GRAININFO_ARR,nodegrainmap,IFsliding_switch)
%
% Computes the element residual contributions due to the two Nitsche terms
% with the averaged stresses and provides an id-array to prepare the
% assembly into the global 'tangentmatrix'.
%
% Input Parameters:
%   xcoords           x-coordinates of current element
%   ycoords           y-coordinates of current element
%   seg_cut_info      information about current subsegment
%   endpoints         endpoints of current interface
%   node              connectivity between all nodes and elements
%   id_dof            mapping between DOFs and enriching grains
%   id_eqns           mappint between nodes and global DOFs (only for this
%                     element)
%   GRAININFO_ARR     some data about the grains
%   nodegrainmap      mapping between nodes and grains
%   IFsliding_switch  indicates sliding case at interface
%
% Returned variables:
%   ke_nit            element stiffnes contribution to global 'tangentmatrix'
%   id                id-array to prepare assembly into 'tangentmatrix'
%

% Author: Matthias Mayr (08/2010)

function [ke_nit id] = lin_nitsche(xcoords,ycoords, ...
  seg_cut_info,endpoints,node,id_dof, ...
  id_eqns,GRAININFO_ARR,nodegrainmap,IFsliding_switch)

%% INITIALIZE
eleID = seg_cut_info.elemno;    % global element ID
elenodes = node(:,eleID);       % global node IDs

intersection = seg_cut_info.xint;

normal = seg_cut_info.normal;                 % normal vector
n_matrix = [normal(1) 0         normal(2);    % normal matrix to keep dimensions 
            0         normal(2) normal(1)];   %   consistent in discrete space
ntn = normal * normal';                       % "normal tensor normal" to project
                                              %   into normal direction
          
pos_g = seg_cut_info.positive_grain;            % positive grain
neg_g = seg_cut_info.negative_grain;            % negative grain
pn_nodes = get_positive_new(eleID,pos_g,neg_g); % pos. and neg. nodes

enrich1 = zeros(1,3);
enrich2 = zeros(1,3);

flg = zeros(1,6);           % flags to evalueat jump in shape functions

ke_nit_12 = zeros(6,12);    % upper right submatrix
ke_nit_21 = zeros(12,6);    % bottom left submatrix
ke_nit_22 = zeros(12,12);   % bottom right submatrix
% ----------------------------------------------------------------------- %
%% GET DISCRETE CONSTITUTIVE MATRICES 'C1' and 'C2'
% grain 1
E = GRAININFO_ARR(pos_g).youngs;
pr = GRAININFO_ARR(pos_g).poisson;
fac = E/(1 - (pr)^2);
C1 = fac * [1.0 pr  0;
            pr  1.0 0.0;
            0   0   (1.-pr)/2 ];

% grain 2
E = GRAININFO_ARR(neg_g).youngs;
pr = GRAININFO_ARR(neg_g).poisson;
fac = E/(1 - (pr)^2);
C2 = fac * [1.0 pr  0;
            pr  1.0 0.0;
            0   0   (1.-pr)/2 ];
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
%% ESTABLISH A SET OF FLAGS
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

% 2 gauss points on segments, since the product of linear shape function
% matrices is integrated, i.e. the integrand is quadratic, so at least 2
% gauss points are required for an exact integration.
gauss = [-sqrt(3)/3 sqrt(3)/3];
weights = [1 1];  % weights can be dropped for two gauss points since
                    % they hace no influence

% get area of element
Area = det([[1 1 1]' xcoords' ycoords'])/2;
% ----------------------------------------------------------------------- %
%% LOOP OVER GAUSS POINTS
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
    N(1,2*c-1:2*c-1)  = N(1,2*c-1:2*c-1)*flg(c);
    N(2,2*c:2*c)      = N(2,2*c:2*c)    *flg(c);
  end;
  
  % compute contributions due to normal constraints
  ke_nit_12 = ke_nit_12 ...
            + CBhat_avg' * n_matrix' * ntn * N * seg_jcob * weights(g);
  ke_nit_21 = ke_nit_21 ...
            + N' * ntn * n_matrix * CBhat_avg * seg_jcob * weights(g);
  ke_nit_22 = ke_nit_22 ...
            + N' * ntn * n_matrix * CBtilde_avg * seg_jcob * weights(g)...
            + CBtilde_avg' * n_matrix' * ntn * N * seg_jcob * weights(g);
          
  % add contributions due to tangential constraints, if interface is fully
  % tied
  if IFsliding_switch == 0
    ke_nit_12 = ke_nit_12 ...
              + CBhat_avg' * n_matrix' * (eye(2) - ntn) * N * seg_jcob * weights(g);
    ke_nit_21 = ke_nit_21 ...
              + N' * (eye(2) - ntn) * n_matrix * CBhat_avg * seg_jcob * weights(g);
    ke_nit_22 = ke_nit_22 ...
              + N' * (eye(2) - ntn) * n_matrix * CBtilde_avg * seg_jcob * weights(g) ...
              + CBtilde_avg' * n_matrix' * (eye(2) - ntn) * N * seg_jcob * weights(g);
  end;
end;
% ----------------------------------------------------------------------- %
%% BUILD ELEMENT STIFFNESS CONTRIBUTION FOR NITSCHE
ke_nit = -1 * [ zeros(6,6)  ke_nit_12;
                ke_nit_21   ke_nit_22];
% if norm(ke_nit - ke_nit',inf) ~= 0
%   ke_nit
%   error('Nitsche-stiffnes not symmetric');
% end;
% ----------------------------------------------------------------------- %
%% BUILD ID-ARRAY FOR ASSEMBLY
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

