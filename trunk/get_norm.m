function [norm] = get_norm(i,j)

global INT_INTERFACE

% Define the z vector out of plane;
z = [0 0 1];

% end points of segment in question
p1 = [INT_INTERFACE(i).shared(2*j-1,1) INT_INTERFACE(i).shared(2*j-1,2)];
p2 = [INT_INTERFACE(i).shared(2*j,1) INT_INTERFACE(i).shared(2*j,2)];

% p is the vector representing the segment
p = [(p2(1) - p1(1)) (p2(2) - p1(2)) 0];
% the normal is the cross product between the two
norm(1) = z(2)*p(3) - z(3)*p(2);
norm(2) = z(3)*p(1) - z(1)*p(3);

%Make norm a unit vector
norm = norm./sqrt(norm(1)^2+norm(2)^2);
