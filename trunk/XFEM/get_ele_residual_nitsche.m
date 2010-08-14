% get_ele_residual_nitsche.m
%
% CALL: get_ele_residual_nitsche()
%
% Computes the element residual contributions due to the two Nitsche terms
% with the averaged stresses.
%
% Input Parameters:
%
%
% Returned variables:
%
%

% Author: Matthias Mayr (08/2010)

function [res_nit id] = get_ele_residual_nitsche(xcoords,ycoords, ...
  seg_cut_info,endpoints,node,x,y,dis,old_ndisp,id_dof,cutlist, ...
  maxngrains,fdisp,id_eqns,GRAININFO_ARR,nodegrainmap)

%% INITIALIZE
eleID = seg_cut_info.elemno;    % global element ID
elenodes = node(:,eleID);       % global node IDs

intersection = seg_cut_info.xint;

normal = seg_cut_info.normal;                 % normal vector
n_matrix = [normal(1) 0         normal(2);    % normal matrix to keep dimensions 
            0         normal(2) normal(1)];   %   consistent in discrete space
          
pos_g = seg_cut_info.positive_grain;            % positive grain
neg_g = seg_cut_info.negative_grain;            % negative grain
pn_nodes = get_positive_new(eleID,pos_g,neg_g); % pos. and neg. nodes

enrich1 = zeros(1,3);
enrich2 = zeros(1,3);

flg = zeros(1,6);

res_nit_base = zeros(6,1);
res_nit_enriched = zeros(12,1);
% ----------------------------------------------------------------------- %
%% GET DISCRETE CONSTITUTIVE MATRICES 'C1' and 'C2'
% grain 1
E = GRAININFO_ARR(pos_g).youngs;
pr = GRAININFO_ARR(pos_g).poisson;
fac = E/(1 - (pr)^2);
C1 = fac*[1.0,  pr,   0;
          pr,   1.0,  0.0;
          0,    0,    (1.-pr)/2 ];

% grain 2
E = GRAININFO_ARR(neg_g).youngs;
pr = GRAININFO_ARR(neg_g).poisson;
fac = E/(1 - (pr)^2);
C2 = fac*[1.0,  pr,   0;
          pr,   1.0,  0.0;
          0,    0,    (1.-pr)/2 ];
        
% clear some temporary variables
clear E pr fac;
% ----------------------------------------------------------------------- %
%% GET AVERAGED STRESS
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

% get stress in current element
[~,stresse] = post_process_better(node,x,y,seg_cut_info.elemno,dis, ...
                old_ndisp,id_dof,cutlist,maxngrains);
              
% stress in positive part of the element (Voigt-notation)
stress1 = [ stresse(1,4,pos_g);
            stresse(1,5,pos_g);
            stresse(1,6,pos_g)];
          
          
% stress in negative part of the element (Voigt-notation)
stress2 = [ stresse(1,4,neg_g);
            stresse(1,5,neg_g);
            stresse(1,6,neg_g)];

% compute averaged stress (Voigt-notation)
stress_avg = 0.5 * (stress1 + stress2);
% ----------------------------------------------------------------------- %
%% GET AVERAGED 'CBhat' and 'CBtilde'
% Since the matrices 'Bhat' and 'Btilde' are constant in an element due to 
% linear shape functions, these quantities can be computed outside the loop
% over the gauss points.

% get B-matrices
[Bhat Btilde1 Btilde2] = getBmatrices(xcoords,ycoords, ...
  elenodes,nodegrainmap,pos_g);


% Since 'Bhat1' = 'Bhat2', <CBhat> = <C>Bhat
CBhat_avg = 0.5 * (C1 + C2) * Bhat;

% Since 'Btilde1' ~= 'Btilde2', the average of 'Btilde' has to be
% considered, too.
CBtilde_avg = 0.5 * (C1 * Btilde1 + C2 * Btilde2);
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
% weights = [1 1];  % weights can be dropped for two gauss points since
                    % they hace no influence

% get area of element
Area = det([[1 1 1]' xcoords' ycoords'])/2;
% ----------------------------------------------------------------------- %
%% LOOP OVER GAUSS POINTS
for g=1:2
  % reset 'N'
  N = zeros(2,12);  % shape function matrix for base DOFs and two possible enrichments
  
  % Get real coordinates of gauss points
  xn = 0.5*(1-gauss(g))*p1(1)+0.5*(1+gauss(g))*p2(1);
  yn = 0.5*(1-gauss(g))*p1(2)+0.5*(1+gauss(g))*p2(2);

  % Evaluate shape functions and assemble 'N'
  for b = 1:3     
    % Get coorindates of area opposite node of concern
    xes = xcoords;
    yes = ycoords;
    xes(b) = xn; 
    yes(b) = yn;
    Larea = det([[1 1 1]' xes' yes'])/2;

    % Evaluate shape function
    N(1,2*b-1)  = N(1,2*b-1)  + Larea/Area;   % First enrichment
    N(2,2*b)    = N(2,2*b)    + Larea/Area;
    N(1,2*b+5)  = N(1,2*b+5)  + Larea/Area;   % Second enrichment
    N(2,2*b+6)  = N(2,2*b+6)  + Larea/Area;
  end;
  
  % build jump by multiplying with 'flg'
  for i=1:6
    N(1,2*i-1)  = N(1,2*i-1)  * flg(i);
    N(2,2*i)    = N(2,2*i)    * flg(i);
  end;
  
  % evaluate gap at current gauss point
  gap = evaluate_gap_gp(xn,yn,flg,xcoords,ycoords,fdisp', ...
    id_eqns,Area);
        
  res_nit_base      = res_nit_base ...
                    + CBhat_avg' * n_matrix' * gap * seg_jcob;
  res_nit_enriched  = res_nit_enriched ...
                    + N' * n_matrix * stress_avg * seg_jcob ...
                    + CBtilde_avg' * n_matrix' * gap * seg_jcob;
end;
% ----------------------------------------------------------------------- %
%% ASSEMBLE ELEMENT RESIDUAL FOR NITSCHE
res_nit = -1 * [res_nit_base;
                res_nit_enriched];
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

% %% INITIALIZE
% % res_nit = zeros(18,1);
% 
% intersection = seg_cut_info.xint;
% 
% enrich1 = zeros(1,3);             % Is there a first enrichment?
% enrich2 = zeros(1,3);             % Is there a second enrichment?
% 
% grain1 = seg_cut_info.grains(1);  % grains, associated to interface
% grain2 = seg_cut_info.grains(2);
% 
% normal = seg_cut_info.normal; % get normal to the interface
% 
% eleID = seg_cut_info.elemno;
% elenodes = node(:,eleID);
% 
% gap = [0 0]';
% 
% B1 = zeros(3,18);
% B2 = zeros(3,18);
% 
% flg = [0 0 0 0 0 0];
% 
% % Establish which nodes are "postively" enriched, and
% % which reside in the "negative" grain
% pos_g = seg_cut_info.positive_grain;
% neg_g = seg_cut_info.negative_grain;
% 
% [pn_nodes] =... 
%    get_positive_new(seg_cut_info.elemno,pos_g,neg_g);
%  
% % res_nit_1 = zeros(1,18);
% % res_nit_2 = zeros(1,18);
% % ----------------------------------------------------------------------- %
% %% GET CONSTITUTIVE MATRICES
% % grain 1
% young = GRAININFO_ARR(pos_g).youngs;
% pr = GRAININFO_ARR(pos_g).poisson;
% fac = young/(1 - (pr)^2);
% C1 = fac*[1.0,  pr,   0;
%           pr,   1.0,  0.0;
%           0,    0,    (1.-pr)/2 ];
% 
% % grain 2
% young = GRAININFO_ARR(neg_g).youngs;
% pr = GRAININFO_ARR(neg_g).poisson;
% fac = young/(1 - (pr)^2);
% C2 = fac*[1.0,  pr,   0;
%           pr,   1.0,  0.0;
%           0,    0,    (1.-pr)/2 ];
% 
% % get cijkl positive (1)
% cijkl_p = find_cijkl(pos_g);
% 
% % get cijkl negative (2)
% cijkl_n = find_cijkl(neg_g);
% % ----------------------------------------------------------------------- %
% %% FIRST PART
%   %% COMPUTE STRESS
%   % Due to linear shape functions, the stress is constant in an element and
%   % so it is not necessary to evaluate it at each gauss point since it is the
%   % same everywhere.
% 
%   % get stress in element
%   [~,stresse] = post_process_better(node,x,y,seg_cut_info.elemno,dis, ...
%     orig_ndisp,id_dof,cutlist,maxngrains);
% 
%   % stress tensors in cut element 'e' on both sides of
%   % interface 'i'
%   sigma1 = [stresse(1,4,grain1) stresse(1,6,grain1);
%             stresse(1,6,grain1) stresse(1,5,grain1)];
%   sigma2 = [stresse(1,4,grain2) stresse(1,6,grain2);
%             stresse(1,6,grain2) stresse(1,5,grain2)];
% 
%   % compute averaged stress in tensor notation
%   sigma_avg = 0.5 * (sigma1 + sigma2);
% 
% %   % put averaged stress in Voigt vector notation
% %   sigma_vec = [ sigma_avg(1,1); 
% %                 sigma_avg(2,2); 
% %                 sigma_avg(1,2)]
% % 
%   % build normal-matrix
%   n_matrix = [normal(1) 0         normal(2);
%               0         normal(2) normal(1)];
% 
%   % compute sigma * n in discrete notation
%   sigmanormal = sigma_avg * normal;%n_matrix * sigma_vec;
%   % --------------------------------------------------------------------- %
%   %% ESTABLISH A SET OF FLAGS
%   % First enrichment
%   for n = 1:3     % loop over nodes
%     % Get the "first" enrichment of node
%     enrich1(n) = id_dof(elenodes(n),3);
% 
%     if enrich1(n) == pos_g
%       if (pn_nodes(n,1) == 1)
%         flg(n) = 1;
%       else
%         flg(n) = 0;
%       end
%     elseif enrich1(n) == neg_g
%       if (pn_nodes(n,2) == 1)
%         flg(n) = -1;
%       else
%         flg(n) = 0;
%       end        
%     end
%   end
% 
%   % Second Enrichment
%   for n = 1:3     % loop over nodes
%     % Get the "second" enrichment of nodes
%     enrich2(n) = id_dof(elenodes(n),5);  
% 
%     if enrich2(n) == pos_g  % If this enrichment corresponds 
%                             % to the positive grain
%       if (pn_nodes(n,1) == 1)
%         flg(3 + n) = 1;
%       else
%         flg(3 + n) = 0;
%       end
% 
%     elseif enrich2(n) == neg_g
%       if (pn_nodes(n,2) == 1)
%         flg(3 + n) = -1;
%       else
%         flg(3 + n) = 0;
%       end        
%     end
%   end
%   % --------------------------------------------------------------------- %
%   %% PREPARE GAUSS QUADRATURE
%   % end points of intersection - direction doesn't matter - this is for the
%   % segment jacobian calculation
% 
%   if all(size(intersection) == [2 2])
%     p1 = intersection(1,:);
%     p2 = intersection(2,:);
%   elseif all(size(intersection) == [1 2])
%     p1 = intersection(1,:);
% 
%     % Second endpoint of segment is also end point of interface
%     endpoint = endpoints(1,:);
% 
%     inside = polygon_contains_point_2d ( 3, [xcoords;ycoords], endpoint );
% 
%     if inside
%       p2 = endpoint;
%     else
%       p2 = endpoints(2,:);
%     end
%   end
% 
%   % jacobian of segment to global
%   he = sqrt((p1(1)-p2(1))^2 + (p1(2)-p2(2))^2);
%   seg_jcob = he/2;
% 
%   % 2 gauss points on segments, since the product of linear shape function
%   % matrices is integrated, i.e. the integrand is quadratic, so at least 2
%   % gauss points are required for an exact integration.
%   gauss = [-sqrt(3)/3 sqrt(3)/3];
% %   weights = [1 1];
% 
%   % get area of element
%   Area = det([[1 1 1]' xcoords' ycoords'])/2;
%   % --------------------------------------------------------------------- %
%   %% LOOP OVER GAUSS POINTS TO ASSEMBLE 'N'
%   % reset shape function matrix 'N'
%   N = zeros(2,18);
%   for g = 1:length(gauss)
%     % Get real coordinates of gauss points
%     xn = 0.5*(1-gauss(g))*p1(1)+0.5*(1+gauss(g))*p2(1);
%     yn = 0.5*(1-gauss(g))*p1(2)+0.5*(1+gauss(g))*p2(2);
% 
%     % Evaluate shape functions and assemble 'N'
%     for b = 1:3     
%       % Get coorindates of area opposite node of concern
%       xes = xcoords;
%       yes = ycoords;
%       xes(b) = xn; 
%       yes(b) = yn;
%       Larea = det([[1 1 1]' xes' yes'])/2;
% 
%       % Evaluate shape function
%       N(1,2*b-1) = N(1,2*b-1) + Larea/Area;   % Base DOFs
%       N(2,2*b)   = N(2,2*b)   + Larea/Area;
%       N(1,2*b+5) = N(1,2*b+5) + Larea/Area;   % First enrichment
%       N(2,2*b+6) = N(2,2*b+6) + Larea/Area;
%       N(1,2*b+11) = N(1,2*b+11) + Larea/Area;   % Second enrichment
%       N(2,2*b+12) = N(2,2*b+12) + Larea/Area;
%     end;
%   end;
% 
%   % Finish gauss quadrature by mutliplying with jacobian and gauss weight
%   % (This can be done outside the loop only because there are only two gauss
%   % points with the same weights = 1).
%   N = N * seg_jcob;
% 
%   % multiply values belonging to nodes in the 'negative' grain with '-1'
%   for c = 1:6
%     N(:,2*c-1+6:2*c+6) = N(:,2*c-1+6:2*c+6)*flg(c);
%   end;
%   % --------------------------------------------------------------------- %
%   %% COMPUTE RESIDUAL PART 1
%   res_nit_1 = -1 * N' * sigmanormal;
%   % --------------------------------------------------------------------- %
% %% SECOND PART
%   %% LOOP OVER GAUSS POINTS TO INTEGRATE OVER GAP
%   for g = 1:2
%     % Get real coordinates of gauss points
%     xn = 0.5*(1-gauss(g))*p1(1)+0.5*(1+gauss(g))*p2(1);
%     yn = 0.5*(1-gauss(g))*p1(2)+0.5*(1+gauss(g))*p2(2);
% 
%     gap = gap + evaluate_gap_penalty_gp(xn,yn,flg,xcoords,ycoords, ...
%                   fdisp',id_eqns,Area);
%   end;
%   % --------------------------------------------------------------------- %
%   %% GET 'B'-MATRICES
%   % Derivatives are constant, so quadrature need not occur over these points
%   % compute derivatives of shape functions in reference coordinates
%   NJr(1) =  1;
%   NJr(2) =  0;
%   NJr(3) = -1;
%   NJs(1) =  0;
%   NJs(2) =  1;
%   NJs(3) = -1;
% 
%   % compute derivatives of x and y wrt psi and eta
%   xr = NJr*xcoords'; % derivatve of x-coord wrt 
%   yr = NJr*ycoords'; 
%   xs = NJs*xcoords';  
%   ys = NJs*ycoords';
% 
%   Jinv = [ys, -yr; -xs, xr];
%   elem_jcob = xr*ys - xs*yr;
% 
%   % compute derivatives of shape functions in element coordinates
%   NJdrs = [NJr; NJs];             % in parameter space
%   NJdxy = Jinv * NJdrs / elem_jcob;   % in real space
% 
%   % Assemble *term 1* (derivatives to go with Cijkl 1)  Positive term
%   NJdxy1 = zeros(2,9);    % derivated shape functions of grain (1)
%   NJdxy1(:,1:3) = NJdxy;
%   for m = 1:3
%     if (pn_nodes(m,1) == 1)  % If the node is enriched positively
%       % Is this the first enrichment?
%       if id_dof(elenodes(m),3) == pos_g
%         % Fill the appropriate slot
%         NJdxy1(:,m+3) = NJdxy(:,m);
%       else  % If this is the second enrichment
%         % Fill the appropriate slot
%         NJdxy1(:,m+6) = NJdxy(:,m);
%       end
%     end
%   end
%   % Now, 'NJdx1' contains the derivatives of the shape functions for grain
%   % (1). First row corresponds to x, second row to y. The zeros in between 
%   % are ommitted here. First 3 colums correspond to base DOFs, second 3 
%   % columns to first enrichment, last 3 columns to a possible second 
%   % enrichment.
% 
%   % Assemble *term 2* (derivitives to go with Cijkl 2)  Negative term
%   NJdxy2 = zeros(2,9);    % derivated shape functions of grain (2)
%   NJdxy2(:,1:3) = NJdxy;
%   for m = 1:3
%     if (pn_nodes(m,2) == 1)  % If the node is enriched negatively
%       % Is this the first enrichment?
%       if id_dof(elenodes(m),3) == neg_g
%         % Fill the appropriate slot
%         NJdxy2(:,m+3) = NJdxy(:,m);
%       else  % If this is the second enrichment
%         % Fill the appropriate slot
%          NJdxy2(:,m+6) = NJdxy(:,m);
%       end
%     end
%   end
%   % Now, 'NJdx2' contains the derivatives of the shape functions for grain
%   % (2). First row corresponds to x, second row to y. The zeros in between 
%   % are ommitted here. First 3 colums correspond to base DOFs, second 3 
%   % columns to first enrichment, last 3 columns to a possible second 
%   % enrichment.
%   
% %   NJdxy1
% %   NJdxy2
%   
%   % Insert zeros in every second slot in 'NJdxy1' and 'NJdxy2'
%   NJ0dxy1 = zeros(2,18);
%   NJ0dxy1 = zeros(2,18);
%   for i=1:9
%     NJ0dxy1(1,2*i-1)  = NJdxy1(1,i);
%     NJ0dxy1(2,2*i)    = NJdxy1(2,i);
%     NJ0dxy2(1,2*i-1)  = NJdxy2(1,i);
%     NJ0dxy2(2,2*i)    = NJdxy2(2,i);
%   end;
%   % --------------------------------------------------------------------- %
%   %% BUILD 'B'-MATRICES
%   % Assing values from 'NJ0dxy1' and 'NJßdxy2' to 'B1' and 'B2' depending
%   % on the grain in which a node resides. Assign value from 'NJ0dxy1' if
%   % the node resides in grain 1 and from 'NJ0dxy2' else.
%   
%   % base DOFs are equivalent
%   % loop over nodes
%   for i=1:3
%     B1(1,2*i-1) = NJ0dxy1(1,2*i-1);
%     B1(2,2*i)   = NJ0dxy1(2,2*i);
%     B1(3,2*i-1) = NJ0dxy1(2,2*i);
%     B1(3,2*i)   = NJ0dxy1(1,2*i-1);
% 
%     B1(1,2*i-1) = NJ0dxy1(1,2*i-1);
%     B1(2,2*i)   = NJ0dxy1(2,2*i);
%     B1(3,2*i-1) = NJ0dxy1(2,2*i);
%     B1(3,2*i)   = NJ0dxy1(1,2*i-1);
%   end;
%   
%   % loop over nodes
%   for i=1:3
%     if nodegrainmap(elenodes(i)) == pos_g
%       % B1
%       B1(1,2*i-1+6) = NJ0dxy1(1,2*i-1+6);
%       B1(2,2*i+6)   = NJ0dxy1(2,2*i+6);
%       B1(3,2*i-1+6) = NJ0dxy1(2,2*i+6);
%       B1(3,2*i+6)   = NJ0dxy1(1,2*i-1+6);
%       
%       B1(1,2*i-1+12) = NJ0dxy1(1,2*i-1+12);
%       B1(2,2*i+12)   = NJ0dxy1(2,2*i+12);
%       B1(3,2*i-1+12) = NJ0dxy1(2,2*i+12);
%       B1(3,2*i+12)   = NJ0dxy1(1,2*i-1+12);
%       
%       % B2
%       B2(1,2*i-1+6) = NJ0dxy2(1,2*i-1+6);
%       B2(2,2*i+6)   = NJ0dxy2(2,2*i+6);
%       B2(3,2*i-1+6) = NJ0dxy2(2,2*i+6);
%       B2(3,2*i+6)   = NJ0dxy2(1,2*i-1+6);
%       
%       B2(1,2*i-1+12) = NJ0dxy2(1,2*i-1+12);
%       B2(2,2*i+12)   = NJ0dxy2(2,2*i+12);
%       B2(3,2*i-1+12) = NJ0dxy2(2,2*i+12);
%       B2(3,2*i+12)   = NJ0dxy2(1,2*i-1+12);
%     else
%       % B1
%       B1(1,2*i-1+6) = NJ0dxy2(1,2*i-1+6);
%       B1(2,2*i+6)   = NJ0dxy2(2,2*i+6);
%       B1(3,2*i-1+6) = NJ0dxy2(2,2*i+6);
%       B1(3,2*i+6)   = NJ0dxy2(1,2*i-1+6);
%       
%       B1(1,2*i-1+12) = NJ0dxy2(1,2*i-1+12);
%       B1(2,2*i+12)   = NJ0dxy2(2,2*i+12);
%       B1(3,2*i-1+12) = NJ0dxy2(2,2*i+12);
%       B1(3,2*i+12)   = NJ0dxy2(1,2*i-1+12);
%       
%       % B2
%       B2(1,2*i-1+6) = NJ0dxy1(1,2*i-1+6);
%       B2(2,2*i+6)   = NJ0dxy1(2,2*i+6);
%       B2(3,2*i-1+6) = NJ0dxy1(2,2*i+6);
%       B2(3,2*i+6)   = NJ0dxy1(1,2*i-1+6);
%       
%       B2(1,2*i-1+12) = NJ0dxy1(1,2*i-1+12);
%       B2(2,2*i+12)   = NJ0dxy1(2,2*i+12);
%       B2(3,2*i-1+12) = NJ0dxy1(2,2*i+12);
%       B2(3,2*i+12)   = NJ0dxy1(1,2*i-1+12);
%     end;
%   end;
%   % --------------------------------------------------------------------- %
%   %% COMPUTE SECOND PART OF RES^{NITSCHE}
%   for g = 1:2
%     % Get real coordinates of gauss points
%     xn = 0.5*(1-gauss(g))*p1(1)+0.5*(1+gauss(g))*p2(1);
%     yn = 0.5*(1-gauss(g))*p1(2)+0.5*(1+gauss(g))*p2(2);
%   
%     gap = gap + evaluate_gap_penalty_gp(xn,yn,flg,xcoords,ycoords, ...
%                     fdisp',id_eqns,Area);
%   end;
%   
%   avgCB = 0.5 *(C1 * B1 + C2 * B2);
%   
%   res_nit_2 = avgCB' * n_matrix' * gap;
% %   for g = 1:2
% %     % Get real coordinates of gauss points
% %     xn = 0.5*(1-gauss(g))*p1(1)+0.5*(1+gauss(g))*p2(1);
% %     yn = 0.5*(1-gauss(g))*p1(2)+0.5*(1+gauss(g))*p2(2);
% %   
% %     gap = evaluate_gap_penalty_gp(xn,yn,flg,xcoords,ycoords, ...
% %                     fdisp',id_eqns,Area);
% %     for i=1:2
% %       for j=1:2
% %         for k=1:2
% %           for l=1:2
% %             for b=1:18
% %               res_nit_2(b) = res_nit_2(b) ...
% %                 + gap(i) * 0.5 * cijkl_p(i,j,k,l) * NJ0dxy1(k,b) * normal(j) ...
% %                 + gap(i) * 0.5 * cijkl_n(i,j,k,l) * NJ0dxy2(k,b) * normal(j);
% %             end;  % b
% %           end;  % l
% %         end;  % k
% %       end;  % j
% %     end;  % i
% %   end;  % g       
%   % --------------------------------------------------------------------- %
% %% ADD BOTH PARTS
% a = res_nit_1'
% b = res_nit_2'
% 
% res_nit = res_nit_1 + res_nit_2;
% end

