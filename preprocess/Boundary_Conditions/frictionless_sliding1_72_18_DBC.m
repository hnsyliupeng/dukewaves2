% frictionless_sliding1_72_18_DBC.m
%
% Dirichlet boundary conditions for frictionless sliding problem based on 
% L = 16, set of structured mesh.  This file is for Dirichlet BCs for 
% 18 elements through the thickness and 72 elements over the length.
%

% Left hand side - displacement boundary conditions

ubar = zeros(2,numnod);
dispbc = zeros(2,numnod);

% x fixed displacements
for i=1:19
    dispbc(1,i) = 1;
%     dispbc(2,i) = 1;
end;

% y fixed displacements
dispbc(2,10) = 1;    % Center on left edge of domain


% Right hand side - displacement boundary conditions
% dispbc(2,447) = 1;  % Center on right edge of domain
dispbc(2,1378) = 1;  % Center on right edge of domain
% dispbc(2,704+19) = 1;
% dispbc(2,704) = 1;
% dispbc(2,704-19) = 1;

% % displacement of right boundary
for i=1369:1387
  dispbc(1,i) = 1;
  ubar(1,i) = -0.4;
end;