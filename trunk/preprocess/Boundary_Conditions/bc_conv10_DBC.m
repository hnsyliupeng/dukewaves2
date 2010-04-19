% bc_conv10_DBC.m
%
% Boundary conditions for beam bending problem based on L = 16, set of
% structured mesh.  This file is for Dirichlet BCs for 50 elements through 
% the thickness and 200 elements over the length.
%

% Left hand side - displacement boundary conditions

% x fixed displacements
dispbc(1,1)  = 1;    % Upper corner
dispbc(1,26) = 1;    % Center
dispbc(1,51) = 1;    % Lower corner

% y fixed displacements
dispbc(2,26) = 1;    % Center