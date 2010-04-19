% frictionless_sliding1_DBC.m
%
% Dirichlet boundary conditions for frictionless sliding problem based on 
% L = 16, set of structed mesh.  This file is for 6 elements through the 
% thickness and 24 elements over the length.
%

% Left hand side - displacement boundary conditions

% x fixed displacements
dispbc(1,1) = 1;    % Upper corner
dispbc(1,2) = 1;    
dispbc(1,3) = 1;    
dispbc(1,4) = 1;    % Center
dispbc(1,5) = 1;    
dispbc(1,6) = 1;    
dispbc(1,7) = 1;    % Lower corner

% y fixed displacements
dispbc(2,4) = 1;    % Center on left edge of domain


% Right hand side - displacement boundary conditions
dispbc(2,172) = 1;  % Center on right edge of domain