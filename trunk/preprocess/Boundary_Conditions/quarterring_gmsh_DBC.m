% quarterring_gmsh_DBC.m
%
% Boundary conditions for gmsh-mesh.
%

% Left hand side - displacement boundary conditions

% x fixed displacements
dispbc(1,1) = 1;
dispbc(1,5) = 1;
dispbc(1,6) = 1;


% y fixed displacements
dispbc(2,1) = 1;
dispbc(2,5) = 1;
dispbc(2,6) = 1;
