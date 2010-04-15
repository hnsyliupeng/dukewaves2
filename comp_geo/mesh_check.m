function [nconn] = mesh_check(conn,x,y)

% Subroutine checks mesh ordering, and corrects.

xe = [];
ye = [];
nconn = zeros(3,1);
    
% get coordinates of element nodes 
for j=1:3
    je = conn(j); xe(j) = x(je); ye(j) = y(je);
end

% Calculate element area
Area = det([[1 1 1]' xe' ye'])/2;
    
% If area is negative, switch nodes 1 and 3
if Area < 0
    nconn(1) = conn(3);
    nconn(2) = conn(2);
    nconn(3) = conn(1);
else
    nconn = conn;
end

