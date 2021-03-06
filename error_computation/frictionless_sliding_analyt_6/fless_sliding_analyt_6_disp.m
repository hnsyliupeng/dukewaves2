% fless_sliding_analyt_2_disp.m
%
% CALL: fless_sliding_analyt_2_disp(xcoord,ycoord)
%
% Compute analytical solution of displacement for example
% 'frictionless_sliding_analyt_2'. Displacement depends on x and y
%
% Input parameters
%   xcoord            x-coordinate
%   ycoord            y-coordinate
%
% Returned variables
%   displacement      vector with x- and y-displacement
%

% Author: Matthias Mayr (05/2010)

function [displacement] = fless_sliding_analyt_2_disp(xcoord,ycoord)

% define material properties
E = 1000.0;         % Young's modulus
nue = 0;            % Poisson's ration

% define geometry
L = 5;

% define load
p = 1;

% define functions for displacements 'ux' and 'uy' (origin of reference
% frame in the bottom left corner of rectangular domain)
ux = @(x,y) -nue*p/2/E/L*x^2 + p/2/E/L*y^2;
uy = @(x,y) -p/E/L*x*y;

% reference for simulation has its origin on the left side, but in the
% middle of the rectangular domain. So, the y-coordinate has to be
% "transformed".
  ycoord = ycoord + 1;

% calculate displacements at (xcoord|ycoord)
dis_x = feval(ux,xcoord,ycoord);
dis_y = feval(uy,xcoord,ycoord);

displacement = [dis_x dis_y];

end