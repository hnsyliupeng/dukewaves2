% getBmatrices.m
%
% CALL: getBmatrices()
%
% Evaluates the B-matrices 'Bhat', 'Btilde1' and 'Btilde2'.
%
% Input Parameters:
%   xcoords         x-coordinates of element's nodes
%   ycoords         y-coordinates of element's nodes
%   elenodes        global IDs of element's nodes
%   nodegrainmap    mapping between nodes and grains
%   pos_g           global ID of "positive" grain
%   id_dof          mapping between nodes and enriching grains
%   
% Returned variables:
%   Bhat
%   Btilde1
%   Btilde2
%

% Author: Matthias Mayr (08/2010)

function [Bhat Btilde1 Btilde2] = getBmatrices(xcoords,ycoords, ...
  elenodes,nodegrainmap,pos_g,id_dof)

% global  NODAL_ENRICH;

%% INITIALIZE
Bhat = zeros(3,6);
% ----------------------------------------------------------------------- %
%% COMPUTE DERIVATIVES OF SHAPE FUNCTIONS
% derivatives of shape functions in reference coordinates
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

Jinv = [ys, -yr; -xs, xr];  % inverse jacobian
jcob = xr*ys - xs*yr;       % det(J) ... determinant of jacobian

% compute derivatives of shape functions in element coordinates
NJdrs = [NJr; NJs];       % in reference space
NJdxy = Jinv*NJdrs/jcob;  % in real space
% ----------------------------------------------------------------------- %
%% ASSEMBLE 'Bhat'
Bhat(1,1:2:5) = NJdxy(1,1:3);  
Bhat(2,2:2:6) = NJdxy(2,1:3);
Bhat(3,1:2:5) = NJdxy(2,1:3);  
Bhat(3,2:2:6) = NJdxy(1,1:3);
% ----------------------------------------------------------------------- %
%% ASSEMBLE 'Btilde1' AND 'Btilde2'
% copy 'Bhat'
Btilde1 = [Bhat Bhat];   % corresponds to positive_grain
Btilde2 = [Bhat Bhat];   % corresponds to negative_grain

% set entries of these nodes to zero, that reside in the other grain
% for i=1:3 % loop over nodes
%   if nodegrainmap(elenodes(i)) == pos_g
%     Btilde2(1:3,2*i-1:2*i) = 0;
%     Btilde2(1:3,2*i-1+6:2*i+6) = 0;
%   else
%     Btilde1(1:3,2*i-1:2*i) = 0;
%     Btilde1(1:3,2*i-1+6:2*i+6) = 0;
%   end;
% end;

% set entries to zero, if there is no second enrichment
for i=1:3 % loop over nodes
  if id_dof(elenodes(i),5) == 0%NODAL_ENRICH(elenodes(i)).cnt - 1 ~= 2
    Btilde1(1:3,6+2*i-1:6+2*i) = 0;
    Btilde2(1:3,6+2*i-1:6+2*i) = 0;
  end;
end;
% ----------------------------------------------------------------------- %
end

