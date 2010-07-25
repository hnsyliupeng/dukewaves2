% fless_sliding_analyt_8_strain.m
%
% CALL: fless_sliding_analyt_8_strain(xcoord,ycoord)
%
% Compute analytical solution of strain for example
% 'frictionless_sliding_analyt_8'. Strain depends on x and y
%
% Input parameters
%   xcoord            x-coordinate
%   ycoord            y-coordinate
%
% Returned variables
%   strain            vector with xx-, yy- and xy-strain
%

% Author: Matthias Mayr (05/2010)

function [strain] = fless_sliding_analyt_8_strain(xcoord,ycoord)

% define material properties
E = 1000.0;         % Young's modulus
nue_1 = 0.3;            % Poisson's ration

% define geometry
H = 4;

% define load
p = 1;

% Elongations


% strains are constant (xcoord|ycoord)
strain_xx = -2.5e-4;
if x_feta < 5
  strain_yy = 7.5e-5;
else
  strain_yy = 0;
end;  
strain_xy = 0;

strain = [strain_xx strain_yy strain_xy];

end