% fless_sliding_analyt_6_traction.m
%
% CALL: fless_sliding_analyt_6_traction(xcoord,ycoord)
%
% Compute analytical solution of tractions for example
% 'frictionless_sliding_analyt_6'. Traction depends on x and y
%
% Input parameters
%   xcoord            x-coordinate
%   ycoord            y-coordinate
%
% Returned variables
%   traction          value of normal traction
%

% Author: Matthias Mayr (05/2010)

function [traction] = fless_sliding_analyt_6_traction(xcoord,ycoord)

% % define material properties
% E = 1000.0;         % Young's modulus
% nue = 0;            % Poisson's ration

% define geometry
L = 16;
H = 1;

% define load
p = 1;

% define functions for displacements 'ux' and 'uy' (origin of reference
% frame in the bottom left corner of rectangular domain)
tr_normal = @(x,y) 0.4 * y;

% calculate displacements at (xcoord|ycoord)
traction = feval(tr_normal,xcoord,ycoord);

end