% bc_conv1_NBC.m
%
% Boundary conditions for beam bending problem based on L = 16, set of
% structed mesh.  This file is Neumann BCs for 2 elements through the 
% thickness and 8 elements over the length.
%

% Right hand side - parabolic distribution of downward shear, P = -1
force(2,25) = -3/16;
force(2,26) = -5/8;
force(2,27) = -3/16;

% Left hand side - force boundary conditions

% y-forces to balance shear on RHS
force(2,1) = 3/16;
force(2,3) = 3/16;

% x-forces to assure linear distribution of stresses

% empty
