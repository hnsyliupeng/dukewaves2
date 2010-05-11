% nit_stiff.m
%
% Computes the Nitsche contribution to the global stiffnes matrix for
% element 'parent'. Also, an id-array for assembly is computed.
%
% Input arguments:
%   node                mapping between elements and their nodes
%   x                   x-coordinates of all nodes
%   y                   y-coordinates of all nodes
%   parent              global element ID of current element
%   id_eqns             mapping between nodes and global DOFs
%   id_dof              shows, if a node is enriched or not
%   pn_nodes            information about "positive" or "negative"
%                       enrichment
%   pos_g               global ID of positively enriched grain
%   neg_g               gloabl ID of negatively enriched grain
%   intersection        coordinates of the intersection points of current
%                       element (points, where the interface cuts the 
%                       element edges)
%   normal              normal vector to the interface
%   IFsliding_switch    indicates, which kind of sliding is chosen
%
% Returned parameters
%   ke_nit              element "stiffness" matrix for Nitsche contribution
%   id                  id-array to enable assembly into global stiffnes
%                       matrix 'bigk'
%

function [ke_nit,id] =...
    nit_stiff(node,x,y,parent,id_eqns,id_dof,pn_nodes,pos_g,neg_g,...
    normal,intersection,endpoints, IFsliding_switch)

% ----------------------------------------------------------------------- %

% INITIALIZE

ke_nit = zeros(18);

xep = [];   % x-coords of nodes of parent element
yep = [];   % y-coords of nodes of parent element
xes = [];
yes = [];
jes = [];
jep = [];   % node IDs of nodes of parent element
cijkl = []; % constitutive tensor

nodes = node(:,parent);

% Get coordinates of parent element
for m=1:3
    jep = node(m,parent); 
    xep(m) = x(jep); 
    yep(m) = y(jep);
end

% ----------------------------------------------------------------------- %

% MATERIAL PROPERTIES

% get cijkl positive (1)
cijkl_p = find_cijkl(pos_g);

% get cijkl negative (2)
cijkl_n = find_cijkl(neg_g);

% ----------------------------------------------------------------------- %

% ASSEMBLE DERIVITAVES

% Derivatives are constant, so quadrature need not occur over these points
% compute derivatives of shape functions in reference coordinates
NJr(1) = 1;
NJr(2) = 0;
NJr(3) = -1;
NJs(1) = 0;
NJs(2) = 1;
NJs(3) = -1;

% compute derivatives of x and y wrt psi and eta
xr = NJr*xep'; % derivatve of x-coord wrt 
yr = NJr*yep'; 
xs = NJs*xep';  
ys = NJs*yep';

Jinv = [ys, -yr; -xs, xr];
elem_jcob = xr*ys - xs*yr;

% compute derivatives of shape functions in element coordinates
NJdrs = [NJr; NJs];             % in parameter space
NJdxy = Jinv*NJdrs/elem_jcob;   % in real space

% Assemble *term 1* (derivitives to go with Cijkl 1)  Positive term
NJdxy1 = zeros(2,9);    % derivated shape functions of grain (1)
NJdxy1(:,1:3) = NJdxy;
for m = 1:3
  if (pn_nodes(m,1) == 1)  % If the node is enriched positively
    % Is this the first enrichment?
    if id_dof(nodes(m),3) == pos_g
      % Fill the appropriate slot
      NJdxy1(:,m+3) = NJdxy(:,m);
    else  % If this is the second enrichment
      % Fill the appropriate slot
      NJdxy1(:,m+6) = NJdxy(:,m);
    end
  end
end
% Now, 'NJdx1' contains the derivatived of the shape functions for grain
% (1). First row corresponds to x, second row to y. The zeros in between 
% are ommitted here. First 3 colums correspond to base DOFs, second 3 
% columns to first enrichment, last 3 columns to a possible second 
% enrichment.

% Assemble *term 2* (derivitives to go with Cijkl 2)  Negative term
NJdxy2 = zeros(2,9);    % derivated shape functions of grain (2)
NJdxy2(:,1:3) = NJdxy;
for m = 1:3
  if (pn_nodes(m,2) == 1)  % If the node is enriched negatively
    % Is this the first enrichment?
    if id_dof(nodes(m),3) == neg_g
      % Fill the appropriate slot
      NJdxy2(:,m+3) = NJdxy(:,m);
    else  % If this is the second enrichment
      % Fill the appropriate slot
       NJdxy2(:,m+6) = NJdxy(:,m);
    end
  end
end
% Now, 'NJdx2' contains the derivatived of the shape functions for grain
% (2). First row corresponds to x, second row to y. The zeros in between 
% are ommitted here. First 3 colums correspond to base DOFs, second 3 
% columns to first enrichment, last 3 columns to a possible second 
% enrichment.
% ----------------------------------------------------------------------- %

% Establish a set of flags
flg = [0 0 0 0 0 0];

% First enrichment
for n = 1:3     % loop over nodes
  % Get the "first" enrichment of node
  enrich1(n) = id_dof(nodes(n),3);
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
  % Get the "second" enrichment of node
  enrich2(n) = id_dof(nodes(n),5);
  if enrich2(n) == pos_g
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
% INTEGRATE SHAPE FUNCTIONS OVER SEGMENT AND ASSEMBLE NITSCHES TERMS

% end points of intersection - direction doesn't matter - this is for the
% segment jacobian calculation only
if all(size(intersection) == [2 2])
  p1 = intersection(1,:);
  p2 = intersection(2,:);
elseif all(size(intersection) == [1 2])
  p1 = intersection(1,:);

  % Second endpoint of segment is also end point of interface
  endpoint = endpoints(1,:);

  % check, whether 'endpoint' is inside the element
  inside = polygon_contains_point_2d ( 3, [xep;yep], endpoint );

  if inside
    p2 = endpoint;
  else
    p2 = endpoints(2,:);
  end
end

% jacobian of segment to global
he = sqrt((p1(1)-p2(1))^2 + (p1(2)-p2(2))^2);
seg_jcob = he/2;

% Gauss points on segments
gauss = [-sqrt(3)/3 sqrt(3)/3];
weights = [1 1];

% loop over Gauss points to assemble N
N = zeros(1,6);   % first three entries correspond to the first enrichment 
                  % of the three nodes, second three entries to the
                  % possible second enrichment.

for g = 1:2
  % Get real coordinates of gauss points
  xn = 0.5*(1-gauss(g))*p1(1)+0.5*(1+gauss(g))*p2(1);
  yn = 0.5*(1-gauss(g))*p1(2)+0.5*(1+gauss(g))*p2(2);

  for b = 1:3     % Evaluate shape functions
    % Get coorindates of area opposite node of concern
    for m=1:3
      jes = node(m,parent); 
      xes(m) = x(jes); 
      yes(m) = y(jes);
    end

    xes(b) = xn; 
    yes(b) = yn;

    Area = det([[1 1 1]' xep' yep'])/2;
    Larea = det([[1 1 1]' xes' yes'])/2;

    % Evaluate shape function
    N(b)   = N(b) + Larea/Area*seg_jcob;
    N(b+3) = N(b+3) + Larea/Area*seg_jcob;
  end
end

% set all entries of not enriched nodes to zero via 'flg'
for b = 1:6
    N(b) = N(b)*flg(b);
end
% Now, 'N' contains the shape function values evaluated at the 2 gauss
% points. The first 3 entries correspond to first enrichment of the 3
% nodes, the second 3 entries to the possible second enrichment of the 3
% nodes.

% ----------------------------------------------------------------------- %
% Fill stiffness matrix twice, once with terms associated with the positive
% grain and once with terms associated with the negative grain.  In each
% case, we are creating a matrix 12x18 in dimension.

ke_nit_sub = zeros(12,18);

% distinguish between different sliding cases
switch IFsliding_switch
  case 0              % no sliding at all (fully tied problem)
    % Computation via nested for-loops, since cijkl is a 4-tensor
    % First, the positive terms.
    for a = 1:6
      for b = 1:9
        for m = 1:2
          for n = 1:2
            p = 2*(a-1) + m;
            q = 2*(b-1) + n;

            for w = 1:2
              for v = 1:2
                ke_nit_sub(p,q) = ke_nit_sub(p,q) +...
                  N(a)*cijkl_p(m,w,n,v)*NJdxy1(v,b)*normal(w)/2;
              end
            end
          end
        end
      end
    end

    % Second, the negative terms.
    for a = 1:6
      for b = 1:9
        for m = 1:2
          for n = 1:2
            p = 2*(a-1) + m;
            q = 2*(b-1) + n;

            for w = 1:2
              for v = 1:2
                ke_nit_sub(p,q) = ke_nit_sub(p,q) +...
                  N(a)*cijkl_n(m,w,n,v)*NJdxy2(v,b)*normal(w)/2;
              end
            end
          end
        end
      end
    end

  case 1              % frictionless sliding
    % Computation via nested for-loops, since cijkl is a 4-tensor.
    % For frictionless sliding, the 'N' and 'cijkl * strain' have to be
    % doted with the normal. 

    % construct a real matrix 'N'
    Nvec = N;
    clear N;
    N = [Nvec(1) 0       Nvec(2) 0       Nvec(3) 0       Nvec(4) 0       Nvec(5) 0       Nvec(6) 0;
         0       Nvec(1) 0       Nvec(2) 0       Nvec(3) 0       Nvec(4) 0       Nvec(5) 0       Nvec(6)];
    
    % dot 'N' with the normal due to frictionless sliding
    N = normal' * N;
           
    % First, the positive terms.
    for i = 1:12
      for b = 1:9
        for m = 1:2
          for n = 1:2
            a = ceil(i/2);
            p = 2*(a-1) + m;
            q = 2*(b-1) + n;

            for w = 1:2
              for v = 1:2
                ke_nit_sub(p,q) = ke_nit_sub(p,q) +...
                  N(i)*normal(m)*cijkl_p(m,w,n,v)*NJdxy1(v,b)*normal(w)/2;
              end;
            end;
          end;
        end;
      end;
    end;

      % Second, the negative terms.
    for i = 1:12
      for b = 1:9
        for m = 1:2
          for n = 1:2
            a = ceil(i/2);
            p = 2*(a-1) + m;
            q = 2*(b-1) + n;

            for w = 1:2
              for v = 1:2
                ke_nit_sub(p,q) = ke_nit_sub(p,q) +...
                  N(i)*normal(m)*cijkl_n(m,w,n,v)*NJdxy2(v,b)*normal(w)/2;
              end
            end
          end
        end
      end
    end

  case 2              % perfect plasticity
    warning('MATLAB:XFEM:main_xfem',...
        'There exists no code for perfect plasticity, yet.')
  case 3              % frictional contact (Coulomb)
    warning('MATLAB:XFEM:main_xfem',...
      'There exists no code for frictional contact (Coulomb), yet.')
  otherwise
    warning('MATLAB:XFEM:main_xfem',...
      'Unvalid slidingID. Choose valid ID or add additional case to switch-case-structure')
end;

ke_nit = [zeros(6,18);
          ke_nit_sub];

% Build id array
nodes = node(:,parent);
id(1)  = id_eqns(nodes(1),1);  % original x dof
id(2)  = id_eqns(nodes(1),2);  % original y dof
id(3)  = id_eqns(nodes(2),1);  % original x dof
id(4)  = id_eqns(nodes(2),2);  % original y dof
id(5)  = id_eqns(nodes(3),1);  % original x dof
id(6)  = id_eqns(nodes(3),2);  % original y dof
id(7)  = id_eqns(nodes(1),3);  % 1st extra x dof
id(8)  = id_eqns(nodes(1),4);  % 1st extra y dof
id(9)  = id_eqns(nodes(2),3);  % 1st extra x dof
id(10) = id_eqns(nodes(2),4);  % 1st extra y dof
id(11) = id_eqns(nodes(3),3);  % 1st extra x dof
id(12) = id_eqns(nodes(3),4);  % 1st extra y dof
id(13) = id_eqns(nodes(1),5);  % 2nd extra x dof
id(14) = id_eqns(nodes(1),6);  % 2nd extra y dof
id(15) = id_eqns(nodes(2),5);  % 2nd extra x dof
id(16) = id_eqns(nodes(2),6);  % 2nd extra y dof
id(17) = id_eqns(nodes(3),5);  % 2nd extra x dof
id(18) = id_eqns(nodes(3),6);  % 2nd extra y dof