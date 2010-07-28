function [ke_sub] = elemstiff_class_subelement(node,x,y,elemi,parentid,id_dof)

% 2d TRI element stiffness routine. The three nodes may have extra degrees
% of freedom. This is a subelement with a parent.

global GRAININFO_ARR CONN
global NODAL_ENRICH X Y
global SUBELEMENT_GRAIN_MAP
dof = zeros(1,3);

% Determine which grain contains this (sub)element
grain = SUBELEMENT_GRAIN_MAP(elemi);

% Determine how many degrees of freedom each parent node has
for i = 1:3
    dof(i) = 2 + 2*(NODAL_ENRICH(node(i,parentid)).cnt - 1);   % the value stored in "cnt is 
end                                                 % one more than the number of
                                                    % enrichments present. 
tot_dof = dof(1) + dof(2) + dof(3);                                            
ke_sub = zeros(tot_dof,tot_dof);
xep = [];   % parent coords
yep = [];
xes = [];   % subelement coords
yes = [];


% First, the standard degrees of freedom

% plane stress D matrix
young = GRAININFO_ARR(grain).youngs;
pr = GRAININFO_ARR(grain).poisson;
fac = young/(1 - (pr)^2);
D = fac*[1.0, pr, 0;
         pr, 1.0, 0.0;
         0, 0, (1.-pr)/2 ];

% % plane strain D matrix
% young = GRAININFO_ARR(grain).youngs;
% pr = GRAININFO_ARR(grain).poisson;
% lam1 = pr*young/((1+pr)*(1-2*pr));
% lam2 = young/(2*(1+pr));
% D = [lam1+2*lam2, lam1, 0;
%      lam1, lam1+2*lam2, 0.0;
%      0, 0, lam2 ];     
     
% get coordinates of PARENT element nodes 
for j=1:3
    jep = node(j,parentid); xep(j) = x(jep); yep(j) = y(jep);
end

% get coordinates of SUB element nodes
for j = 1:3
    jes = CONN(j,elemi); xes(j) = X(jes); yes(j) = Y(jes);
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
xr = NJr*xep'; yr = NJr*yep'; xs = NJs*xep';  ys = NJs*yep';
Jinv = [ys, -yr; -xs, xr];
jcob = abs(xr*ys - xs*yr);
% compute derivatives of shape functions of PARENT in element coordinates
NJdrs = [NJr; NJs];
NJdxy = Jinv*NJdrs/jcob;
% assemble B matrix
BJ = zeros(3,6);
BJ(1,1:2:5) = NJdxy(1,1:3);  BJ(2,2:2:6) = NJdxy(2,1:3);
BJ(3,1:2:5) = NJdxy(2,1:3);  BJ(3,2:2:6) = NJdxy(1,1:3);
% Area of the SUBelement
Area = abs(det([[1 1 1]' xes' yes'])/2);
% Area = 0.5;
% assemble ke
ke_sub(1:6,1:6) = ke_sub(1:6,1:6) + Area*BJ'*D*BJ;

% Assemble B matrix for enriched nodes

BA = [];
% for each node
for i = 1:3
    extra = (dof(i) - 2)/2;
    Ba = zeros(3,2,extra);
    % for each enrichment
    for j = 1:extra
        % Heaviside, does the enrichment coorespond to the current domain?
        if id_dof(node(i,parentid),(2 + (2*j-1))) == grain;
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

% Assemble element stiffness
ke_sub(1:6,7:tot_dof) = Area*BJ'*D*BA;
ke_sub(7:tot_dof,1:6) = Area*BA'*D*BJ;
ke_sub(7:tot_dof,7:tot_dof) = Area*BA'*D*BA;

