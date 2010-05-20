% fless_sliding_analyt_8_disp.m
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

function [displacement] = fless_sliding_analyt_8_disp(xcoord,ycoord)

% define material properties
E = 1000.0;         % Young's modulus
nue = 0;            % Poisson's ration

% define geometry
H = 4;

% define load
p = 1;

% define functions for displacements 'ux' and 'uy' (origin of reference
% frame in the bottom left corner of rectangular domain)
ux = @(x,y) 2*p/H/E*(x-8)*y;
uy = @(x,y) -nue/E*p/H*y^2 - p/H/E*(x-8)^2;

% calculate displacements at (xcoord|ycoord)
dis_x = feval(ux,xcoord,ycoord);
dis_y = feval(uy,xcoord,ycoord);

displacement = [dis_x dis_y];

end