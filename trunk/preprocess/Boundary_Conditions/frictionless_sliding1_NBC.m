% frictionless_sliding1_NBC.m
%
% Neumann boundary conditions for frictionless sliding problem based on 
% L = 16, set of structed mesh.  This file is for 6 elements through the 
% thickness and 24 elements over the length.
%

% Author: Matthias Mayr (04/2010)

% maximum load in x-direction
maxload = 1000;

% nodal load
nodeload = maxload/6;

% Right hand side - parabolic distribution of downward shear, P = -1
force(1,169) = nodeload/2;
force(1,170) = nodeload;
force(1,171) = nodeload;
force(1,172) = nodeload;
force(1,173) = nodeload;
force(1,174) = nodeload;
force(1,175) = nodeload/2;