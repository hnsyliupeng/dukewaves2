% bc_conv1_NBC.m
%
% Boundary conditions for beam bending problem based on L = 16, set of
% structed mesh.  This file is Neumann BCs for 2 elements through the 
% thickness and 8 elements over the length.
%

% load = 1.0e-12;
load = 100;

% Right hand side - parabolic distribution of downward shear, P = -1
% force(1,7) = load/2;
% force(1,8) = load/2;

force(1,10) = load/4;
force(1,11) = load/2;
force(1,12) = load/4;

