% get_ele_residual_int.m
%
% CALL: get_ele_residual_int(xcoords,ycoords,dis, ...
%                              totaldis,cut)
%
% This function computes the contribution of the internal forces of an
% element due to its deformation (e.g. elastic deformation) to the residual
% of the element
%
% Input parameters:
%   elenodes        global node IDs of current element
%   xcoords         x-coordinates of element's nodes
%   ycoords         y-coordinates of element's nodes
%   dis             global displacement vector (re-assembled)
%   totaldis        global solution vector for element's DOFs
%   cut             entry in 'cutlist' of current element
%
% Returned variables:
%   residual        internal force contribution to the element's residual

% Author: Matthias Mayr (07/2010)

function [residual] = get_ele_residual_int(elenodes,xcoords,ycoords,dis, ...
  orig_ndisp,cut,maxngrains,e,id_dof)

global NODAL_ENRICH
global SUBELEMENT_GRAIN_MAP

% initialize
dof = zeros(1,3);
residual = [];

% Determine which grain contains this element
grain = SUBELEMENT_GRAIN_MAP(e);

% Determine how many degrees of freedom each node has
for i = 1:3
  dof(i) = 2 + 2*(NODAL_ENRICH(elenodes(i)).cnt - 1);   % the value stored in "cnt is 
end                                                   % one more than the number of
                                                      % enrichments present. 
tot_dof = dof(1) + dof(2) + dof(3);                                            

% First, the standard degrees of freedom
% compute derivatives of shape functions in reference coordinates
NJr(1) = 1;
NJr(2) = 0;
NJr(3) = -1;
NJs(1) = 0;
NJs(2) = 1;
NJs(3) = -1;

% compute derivatives of x and y wrt psi and eta
xr = NJr*xcoords'; 
yr = NJr*ycoords'; 
xs = NJs*xcoords';  
ys = NJs*ycoords';
Jinv = [ys, -yr; -xs, xr];
jcob = xr*ys - xs*yr;

% compute derivatives of shape functions in element coordinates
NJdrs = [NJr; NJs];       % in parameter coordinates 'r' and 's'
NJdxy = Jinv*NJdrs/jcob;  % in real coordinates 'x'and 'y'

% assemble B matrix: B = L * N
BJ = zeros(3,6);
BJ(1,1:2:5) = NJdxy(1,1:3);  
BJ(2,2:2:6) = NJdxy(2,1:3);
BJ(3,1:2:5) = NJdxy(2,1:3);  
BJ(3,2:2:6) = NJdxy(1,1:3);

% Area of the element
Area = det([[1 1 1]' xcoords' ycoords'])/2;

% id = [id_eqns(elenodes(1),1:2) id_eqns(elenodes(2),1:2) id_eqns(elenodes(3),1:2)...
%     id_eqns(elenodes(1),3:6) id_eqns(elenodes(2),3:6) id_eqns(elenodes(3),3:6)];
% 
% eliminate = find(id == 0);
% for i = size(eliminate,2):-1:1  % eliminate dofs with index '0'
%     id(eliminate(i)) = [];
% end

% Assemble B matrix for enriched nodes
BA = [];

% if no enrichment, end there
if tot_dof ~= 6

  % for each node
  for i=1:3
    extra = (dof(i) - 2)/2;
    Ba = zeros(3,2,extra);
    if extra ~= 0
      % for each enrichment
      for j = 1:extra
        % Heaviside, does the enrichment coorespond to the current domain?
        if id_dof(i,(2 + (2*j-1))) == grain;
          H = 1;
        else
          H = 0;
        end
        Ba(1,1,j) = NJdxy(1,i) * H;
        Ba(2,2,j) = NJdxy(2,i) * H;
        Ba(3,1,j) = NJdxy(2,i) * H;
        Ba(3,2,j) = NJdxy(1,i) * H;
        BA = [BA Ba(:,:,j)];
      end
    end
  end  
  e,BA
end;

% get stresses in element
[straine stresse] = get_ele_stress(elenodes,xcoords,ycoords,dis,...
  orig_ndisp,BJ,cut,maxngrains,e,id_dof);

residual = zeros(tot_dof,1);

for g = 1:size(stresse,3)
  residual = residual + [BJ BA]' * stresse(1,4:6,g)' * Area;   % base DOFs
end;


end

