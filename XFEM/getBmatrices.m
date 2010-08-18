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
  elenodes,nodegrainmap,pos_g,neg_g,pn_nodes,id_dof)

% global  NODAL_ENRICH;

%% INITIALIZE
Bhat    = zeros(3,6);
Btilde1 = zeros(3,12);
Btilde2 = zeros(3,12);
% ----------------------------------------------------------------------- %
%% COMPUTE DERIVATIVES OF SHAPE FUNCTIONS
% Derivatives are constant, so quadrature need not occur over these points
% compute derivatives of shape functions in reference coordinates
NJr(1) = 1;
NJr(2) = 0;
NJr(3) = -1;
NJs(1) = 0;
NJs(2) = 1;
NJs(3) = -1;

% compute derivatives of x and y wrt psi and eta
xr = NJr*xcoords'; % derivatve of x-coord wrt r
yr = NJr*ycoords'; % derivatve of y-coord wrt r
xs = NJs*xcoords'; % derivatve of x-coord wrt s
ys = NJs*ycoords'; % derivatve of y-coord wrt s

Jinv = [ys, -yr; -xs, xr];  % inverse jacobian
elem_jcob = xr*ys - xs*yr;  % determinant of jacobian

% compute derivatives of shape functions in element coordinates
NJdrs = [NJr; NJs];             % in parameter space
NJdxy = Jinv*NJdrs/elem_jcob;   % in real space

% Assemble *term 1* (derivatives to go with Cijkl 1)  Positive term
NJdxy1 = zeros(2,9);    % derivated shape functions of grain (1)
NJdxy1(:,1:3) = NJdxy;
for m = 1:3
  if (pn_nodes(m,1) == 1)  % If the node is enriched positively
    % Is this the first enrichment?
    if id_dof(elenodes(m),3) == pos_g
      % Fill the appropriate slot
      NJdxy1(:,m+3) = NJdxy(:,m);
    else  % If this is the second enrichment
      % Fill the appropriate slot
      NJdxy1(:,m+6) = NJdxy(:,m);
    end
  end
end
% Now, 'NJdx1' contains the derivatives of the shape functions for grain
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
    if id_dof(elenodes(m),3) == neg_g
      % Fill the appropriate slot
      NJdxy2(:,m+3) = NJdxy(:,m);
    else  % If this is the second enrichment
      % Fill the appropriate slot
       NJdxy2(:,m+6) = NJdxy(:,m);
    end
  end
end
% Now, 'NJdx2' contains the derivatives of the shape functions for grain
% (2). First row corresponds to x, second row to y. The zeros in between 
% are ommitted here. First 3 colums correspond to base DOFs, second 3 
% columns to first enrichment, last 3 columns to a possible second 
% enrichment.
% ----------------------------------------------------------------------- %
%% ASSEMBLE 'Bhat'
% loop over nodes
for i=1:3
  Bhat(1,2*i-1) = NJdxy(1,i);
  Bhat(2,2*i)   = NJdxy(2,i);
  Bhat(3,2*i-1) = NJdxy(2,i);
  Bhat(3,2*i)   = NJdxy(1,i);
end;
% ----------------------------------------------------------------------- %
%% ASSEMBLE 'Btilde1' AND 'Btilde2'
% loop over nodes
for i=1:3
  % first enrichment in positive grain
  Btilde1(1,2*i-1)      = NJdxy1(1,i+3);
  Btilde1(2,2*i)        = NJdxy1(2,i+3);
  Btilde1(3,2*i-1)      = NJdxy1(2,i+3);
  Btilde1(3,2*i)        = NJdxy1(1,i+3);
  
  % second enrichment in positive grain
  Btilde1(1,2*i-1 + 6)  = NJdxy1(1,i+6);
  Btilde1(2,2*i + 6)    = NJdxy1(2,i+6);
  Btilde1(3,2*i-1 + 6)  = NJdxy1(2,i+6);
  Btilde1(3,2*i + 6)    = NJdxy1(1,i+6);
  
  % first enrichment in negative grain
  Btilde2(1,2*i-1)      = NJdxy2(1,i+3);
  Btilde2(2,2*i)        = NJdxy2(2,i+3);
  Btilde2(3,2*i-1)      = NJdxy2(2,i+3);
  Btilde2(3,2*i)        = NJdxy2(1,i+3);
  
  % second enrichment in negative grain
  Btilde2(1,2*i-1 + 6)  = NJdxy2(1,i+6);
  Btilde2(2,2*i + 6)    = NJdxy2(2,i+6);
  Btilde2(3,2*i-1 + 6)  = NJdxy2(2,i+6);
  Btilde2(3,2*i + 6)    = NJdxy2(1,i+6);
end;

% if elenodes == [425 426 445]'
%   pos_g
%   neg_g
%   NJdxy1
%   Btilde1
%   NJdxy2
%   Btilde2
% end;
% ----------------------------------------------------------------------- %
end

