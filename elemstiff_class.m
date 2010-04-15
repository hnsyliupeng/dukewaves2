function [id,ke] = elemstiff_class(node,x,y,e,id_dof,id_eqns)

% 2d TRI element stiffness routine. The three nodes may have extra degrees
% of freedom, but this element is UNCUT.

global GRAININFO_ARR
global NODAL_ENRICH
global SUBELEMENT_GRAIN_MAP
dof = zeros(1,3);

% Determine which grain contains this element
grain = SUBELEMENT_GRAIN_MAP(e);

% Determine how many degrees of freedom each node has
for i = 1:3
    dof(i) = 2 + 2*(NODAL_ENRICH(node(i,e)).cnt - 1);   % the value stored in "cnt is 
end                                                 % one more than the number of
                                                    % enrichments present. 
tot_dof = dof(1) + dof(2) + dof(3);                                            
ke = zeros(tot_dof,tot_dof);
xe = [];
ye = [];

% First, the standard degrees of freedom

% plane stress D matrix
young = GRAININFO_ARR(grain).youngs;
pr = GRAININFO_ARR(grain).poisson;
fac = young/(1 - (pr)^2);
D = fac*[1.0, pr, 0;
         pr, 1.0, 0.0;
         0, 0, (1.-pr)/2 ];
      
% get coordinates of element nodes 
for j=1:3
    je = node(j,e); xe(j) = x(je); ye(j) = y(je);
end

% compute element stiffness
% compute derivatives of shape functions in reference coordinates
NJr(1) = 1;
NJr(2) = 0;
NJr(3) = -1;
NJs(1) = 0;
NJs(2) = 1;
NJs(3) = -1;
% compute derivatives of x and y wrt psi and eta
xr = NJr*xe'; yr = NJr*ye'; xs = NJs*xe';  ys = NJs*ye';
Jinv = [ys, -yr; -xs, xr];
jcob = xr*ys - xs*yr;
% compute derivatives of shape functions in element coordinates
NJdrs = [NJr; NJs];
NJdxy = Jinv*NJdrs/jcob;
% assemble B matrix
BJ = zeros(3,6);
BJ(1,1:2:5) = NJdxy(1,1:3);  BJ(2,2:2:6) = NJdxy(2,1:3);
BJ(3,1:2:5) = NJdxy(2,1:3);  BJ(3,2:2:6) = NJdxy(1,1:3);
% Area of the element
Area = det([[1 1 1]' xe' ye'])/2;
% Area = 0.5;
% assemble ke
ke(1:6,1:6) = ke(1:6,1:6) + Area*BJ'*D*BJ;

id = [id_eqns(node(1,e),1:2) id_eqns(node(2,e),1:2) id_eqns(node(3,e),1:2)...
    id_eqns(node(1,e),3:6) id_eqns(node(2,e),3:6) id_eqns(node(3,e),3:6)];
eliminate = find(id == 0);
for i = size(eliminate,2):-1:1
    id(eliminate(i)) = [];
end

% if no enrichment, end there

if tot_dof == 6
    return
end

% Assemble B matrix for enriched nodes

BA = [];

% for each node
i = 1;
while i < 4
    extra = (dof(i) - 2)/2;
    Ba = zeros(3,2,extra);
    if extra ~= 0
        % for each enrichment
        for j = 1:extra
            % Heaviside, does the enrichment coorespond to the current domain?
            if id_dof(node(i,e),(2 + (2*j-1))) == grain;
                H = 1;
            else
                H = 0;
            end
            Ba(1,1,j) = NJdxy(1,i)*H;
            Ba(2,2,j) = NJdxy(2,i)*H;
            Ba(3,1,j) = NJdxy(2,i)*H;
            Ba(3,2,j) = NJdxy(1,i)*H;
            BA = [BA Ba(:,:,j)];
        end
    end
    i = i + 1;
end  

% Assemble element stiffness
ke(1:6,7:tot_dof) = Area*BJ'*D*BA;
ke(7:tot_dof,1:6) = Area*BA'*D*BJ;
ke(7:tot_dof,7:tot_dof) = Area*BA'*D*BA;
