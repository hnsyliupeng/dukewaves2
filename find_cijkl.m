function [cijkl] = find_cijkl(grn)

global GRAININFO_ARR

% Find Cijkl - a fourth order tensor

grain = grn;

% Grain material properties and D matrix

% plane stress D matrix
young = GRAININFO_ARR(grain).youngs;
pr = GRAININFO_ARR(grain).poisson;
fac = young/(1 - (pr)^2);
D = fac*[1.0, pr, 0;
         pr, 1.0, 0.0;
         0, 0, (1.-pr)/2 ];
     
% Fill Cijkl (4th order tensor)

cijkl(1,1,1,1) = D(1,1);
cijkl(1,1,2,2) = D(1,2);
cijkl(2,2,1,1) = D(2,1);
cijkl(1,1,1,2) = D(1,3);
cijkl(1,1,2,1) = D(1,3);
cijkl(1,2,1,1) = D(3,1);
cijkl(2,1,1,1) = D(3,1);
cijkl(2,2,1,2) = D(2,3);
cijkl(2,2,2,1) = D(2,3);
cijkl(1,2,2,2) = D(3,2);
cijkl(2,1,2,2) = D(3,2);
cijkl(1,2,1,2) = D(3,3);
cijkl(1,2,2,1) = D(3,3);
cijkl(2,1,1,2) = D(3,3);
cijkl(2,1,2,1) = D(3,3);
cijkl(2,2,2,2) = D(2,2);
