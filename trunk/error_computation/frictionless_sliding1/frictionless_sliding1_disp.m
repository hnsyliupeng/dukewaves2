% frictionless_slliding1_disp.m
%
% CALL: fless_sliding_analyt_8_disp(xcoord,ycoord)
%
% Compute analytical solution of displacement for example
% 'frictionless_sliding_analyt_8'. Displacement depends on x and y
%
% Input parameters
%   xcoord            x-coordinate
%   ycoord            y-coordinate
%
% Returned variables
%   displacement      vector with x- and y-displacement
%

% Author: Matthias Mayr (05/2010)

function [displacement] = frictionless_sliding1_disp(xcoord,ycoord)

% define material properties
E = 1000.0;         % Young's modulus
nue1 = 0.3;            % Poisson's ration
nue2 = 0.0;            % Poisson's ration

% define geometry
H = 4;
L = 16;
L1 = 5;
L2 = 11;

% define load
p = -0.25;

% max. displacements of boundaries
deltaL = -4e-3;
deltaH = 3e-4;

% define functions for displacements 'ux' and 'uy' (origin of reference
% frame in the middle of the left boundary of rectangular domain)
ux = @(x,y) deltaL / L * x;
uy = @(x,y) deltaH / H * y;

% calculate displacements at (xcoord|ycoord)
dis_x = feval(ux,xcoord,ycoord);
if xcoord < L1
  dis_y = feval(uy,xcoord,ycoord);
else
  dis_y = 0;
end;

displacement = [dis_x dis_y];

end