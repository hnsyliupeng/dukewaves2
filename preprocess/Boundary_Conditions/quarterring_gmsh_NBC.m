% quarterring_gmsh_NBC.m
%
% Boundary conditions for gmsh-mesh.
%

load = -100;

% Right hand side 
force(2,2) = load;
force(2,3) = load;
force(2,7) = load;

