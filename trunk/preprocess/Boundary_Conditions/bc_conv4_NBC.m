% bc_conv4_NBC.m
%
% Neumann boundary conditions for beam bending problem based on L = 16, set of
% structed mesh.  This file is for 6 elements through the thickness and 24
% elements over the length.
%

% Right hand side - parabolic distribution of downward shear, P = -1
force(2,169) = -11/432;
force(2,170) = -29/216;
force(2,171) = -47/216;
force(2,172) = -53/216;
force(2,173) = -47/216;
force(2,174) = -29/216;
force(2,175) = -11/432;

% Left hand side - force boundary conditions

% y-forces to balance shear on RHS
force(2,1) = 11/432;
force(2,2) = 29/216;
force(2,3) = 47/216;
force(2,5) = 47/216;
force(2,6) = 29/216;
force(2,7) = 11/432;

% x-forces to assure linear distribution of stresses
force(1,2) = -8/3;
force(1,3) = -4/3;
force(1,5) =  4/3;
force(1,6) =  8/3;
