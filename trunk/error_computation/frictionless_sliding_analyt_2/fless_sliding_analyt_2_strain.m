% fless_sliding_analyt_2_strain.m
%
% CALL: fless_sliding_analyt_2_strain(xcoord,ycoord)
%
% Compute analytical solution of strain for example
% 'frictionless_sliding_analyt_2'. Strain depends on x and y
%
% Input parameters
%   xcoord            x-coordinate
%   ycoord            y-coordinate
%
% Returned variables
%   strain            vector with xx-, yy- and xy-strain
%

% Author: Matthias Mayr (05/2010)

function [strain] = fless_sliding_analyt_2_strain(xcoord,ycoord)

% define material properties
E = 1000.0;         % Young's modulus
nue = 0;            % Poisson's ration

% define geometry
L = 5;

% define load
p = 1;

% define functions for displacements 'ux' and 'uy' (origin of reference
% frame in the bottom left corner of rectangular domain)
exx = @(x,y) -nue*p/E/L*x;
eyy = @(x,y) -p/E/L*x;
exy = @(x,y) 0;

% reference for simulation has its origin on the left side, but in the
% middle of the rectangular domain. So, the y-coordinate has to be
% "transformed".
  ycoord = ycoord + 1;

% calculate strains at (xcoord|ycoord)
strain_xx = feval(exx,xcoord,ycoord);
strain_yy = feval(eyy,xcoord,ycoord);
strain_xy = feval(exy,xcoord,ycoord);

strain = [strain_xx strain_yy strain_xy];

end