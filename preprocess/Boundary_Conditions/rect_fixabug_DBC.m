% bc_conv1_DBC.m
%
% Boundary conditions for beam bending problem based on L = 16, set of
% structed mesh.  This file is for the Dirichlet BCs for 2 elements through
% the thickness and 8 elements over the length.
%

% Left hand side - displacement boundary conditions

% x fixed displacements
% dispbc(1,1) = 1;
% dispbc(1,2) = 1;


% y fixed displacements
% dispbc(2,1) = 1;
dispbc(2,2) = 1;


dispbc(1,1) = 1;
dispbc(1,2) = 1;
dispbc(1,3) = 1;
