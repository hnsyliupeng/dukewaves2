%% gen_penalty_alternative.m
%
% CALL: gen_penalty_alternative(node,x,y,parent,id_eqns,id_dof, ...
%         pn_nodes,pos_g,neg_g,intersection,endpoints,normal, ...
%         IFsliding_switch,slidestate,contactstate)
%
% Computes the penalty-contribution to the global stiffnes matrix for
% element 'parent'. Also, an id-array for assembly is computed. This method
% uses the classical approach where only one integral is used to represent
% the penalty contribution. The name "alternative" leans to 'Sanders2009'
% and is only used, since the other method was implemented first.
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
%   slidestate          sliding state of current element (stick or slip)
%   contactstate        current contact state (open = 0, closed = 1)
%
% Returned variables
%   ke_pen              element "stiffness" matrix for penalty contribution
%   id                  id-array to enable assembly into global stiffnes
%                       matrix 'bigk'
%

function [ke_pen,id] =... 
    gen_penalty_alternative(node,x,y,parent,id_eqns,id_dof,...
                 pn_nodes,pos_g,neg_g,intersection,endpoints,normal, ...
                 IFsliding_switch,slidestate,contactstate)

%% Initialize
% penalty stiffness contribution
ke_pen = zeros(12);

% some coordinate variables
xep = [];
yep = [];
xes = [];
yes = [];

% global node IDs of current element
nodes = node(:,parent);

% shape function matrix for first and second enrichments
N = zeros(2,12);

% Get coordinates of parent element
for m=1:3
  jep = node(m,parent); 
  xep(m) = x(jep); 
  yep(m) = y(jep);
end
% ----------------------------------------------------------------------- %
%% Establish a set of flags
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
  % Get the "second" enrichment of nodes
  enrich2(n) = id_dof(nodes(n),5);  

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

% 3 gauss points on segments, since the product of linear shape function
% matrices is integrated, i.e. the integrand is quadratic, so at least 3
% gauss points are required for an exact integration.
gauss = [-sqrt(3/5) 0 sqrt(3/5)];
weights = [5/9 8/9 5/9];
% ----------------------------------------------------------------------- %
%% loop over Gauss points to assemble ke_pen
for g = 1:3
  ke_pen_additional = zeros(12,12);
  
  %% assemble N
  % Get real coordinates of gauss points
  xn = 0.5*(1-gauss(g))*p1(1)+0.5*(1+gauss(g))*p2(1);
  yn = 0.5*(1-gauss(g))*p1(2)+0.5*(1+gauss(g))*p2(2);

  for b = 1:3     % Evaluate shape functions
    % Get coorindates of area opposite node of concern
    for m=1:3
      jes = node(m,parent); 
      xes(m) = x(jes); 
      yes(m) = y(jes);
    end;

    xes(b) = xn; 
    yes(b) = yn;

    Area = det([[1 1 1]' xep' yep'])/2;
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
  
  %% treatment of different sliding cases
  switch IFsliding_switch
    case 0              % no sliding at all (fully constrained)
      % compute 'ke_pen'
      ke_pen_additional = (N'*N);
    case 1              % frictionless sliding
      % dot 'N' with the normal vector 'normal'
      N = N' * normal; 

      % compute ke_pen
      ke_pen_additional = (N * N');
    case 2              % perfect plasticity
      % compute 'ke_pen' depending on current slide state (stick or slip)
      if slidestate == 0      % stick
        % equates to fully tied case
        % compute 'ke_pen'
        ke_pen_additional = (N'*N);    
      elseif slidestate == 1  % slip
        % equates to frictionless sliding case
        % dot 'N' with the normal vector 'normal' 
        N = N' * normal;     

        % compute 'ke_pen'
        ke_pen_addtional = (N * N');  
      else
        error('MATLAB:XFEM:UnvalidSlideState','Unvalid slidestate');
      end;
    case 3              % frictional contact (Coulomb)
      warning('MATLAB:XFEM:main_xfem',...
        'There exists no code for frictional contact (Coulomb), yet.')
    case 4              % frictionless contact (only opening contact)
       % compute 'ke_pen' depending on current slide state (stick or slip)
      if contactstate == 0      % open contact
        % no constraints required
        ke_pen_additional = zeros(size(ke_pen));   
      else                      % closed frictionless contact
        % dot 'N' with the normal vector 'normal'
        N = N' * normal; 

        % compute ke_pen
        ke_pen_additional = (N * N');
      end;
    otherwise
      warning('MATLAB:XFEM:main_xfem',...
        'Unvalid slidingID. Choose valid ID or add additional case to switch-case-structure')
  end;
  
  % finish the gauss quadrature by multiplying with gauss weight and
  % jacobian
  ke_pen = ke_pen + ke_pen_additional .* weights(g) * seg_jcob;
  
  % reinitialize shape function matrix 'N' for next gauss point
  N = zeros(2,12);
  % --------------------------------------------------------------------- %
end;
% ----------------------------------------------------------------------- %
%% Build id array
% get nodes of element
nodes = node(:,parent);

% get DOFs for first enrichment
id(1) = id_eqns(nodes(1),3);  % 1st extra x dof
id(2) = id_eqns(nodes(1),4);  % 1st extra y dof
id(3) = id_eqns(nodes(2),3);  % 1st extra x dof
id(4) = id_eqns(nodes(2),4);  % 1st extra y dof
id(5) = id_eqns(nodes(3),3);  % 1st extra x dof
id(6) = id_eqns(nodes(3),4);  % 1st extra y dof

% get DOFs for second enrichment
id(7)  = id_eqns(nodes(1),5);  % 2nd extra x dof
id(8)  = id_eqns(nodes(1),6);  % 2nd extra y dof
id(9)  = id_eqns(nodes(2),5);  % 2nd extra x dof
id(10) = id_eqns(nodes(2),6);  % 2nd extra y dof
id(11) = id_eqns(nodes(3),5);  % 2nd extra x dof
id(12) = id_eqns(nodes(3),6);  % 2nd extra y dof
% ----------------------------------------------------------------------- %