% frictionless_sliding1_72_18_NBC.m
%
% Neumann boundary conditions for frictionless sliding problem based on 
% L = 16, set of structured mesh.  This file is for Neumann BCs for 
% 18 elements through the thickness and 72 elements over the length.
%


% Author: Matthias Mayr (04/2010)


% maximum load in x-direction
totalload = 0;%-1000;   % 'totalload' = 'constant traction' * 'length of load'


% nodal load
nodeload = totalload/18;


force = zeros(2,numnod);


% Right hand side - constant traction
force(1,1369)=nodeload/2;
for i=1370:1386
    force(1,i)=nodeload;
end;
force(1,1387)=nodeload/2;