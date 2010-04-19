% bc_conv4_DBC.m
%
% Dirichlet boundary conditions for beam bending problem based on L = 16, set of
% structed mesh.  This file is for 6 elements through the thickness and 24
% elements over the length.
%

% Left hand side - displacement boundary conditions

% x fixed displacements
dispbc(1,1) = 1;    % Upper corner
dispbc(1,4) = 1;    % Center
dispbc(1,7) = 1;    % Lower corner

% y fixed displacements
dispbc(2,4) = 1;    % Center