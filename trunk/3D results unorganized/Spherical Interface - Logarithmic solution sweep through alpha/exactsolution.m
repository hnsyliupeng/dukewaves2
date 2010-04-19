function [uex] = exactsolution(x,y,z)

c1 = 1;
c2 = 1;

%%%%% Patch Test Planar Interface %%%%%%%%
% uex = c1 + c2*z;

%%%%%% Spherical Interface %%%%%%%%%%%
rad =  sqrt((x)*(x) + (y)*(y)+ (z)*(z));
% if(rad~=0)
%     uex = c1 + c2*log(rad);
% else
%     uex = 0;
% end

%%%%%%%% Quadratic Solution %%%%%%%%%%%
% uex = (c1 + c2*z*z*z);
%%%%%%%% Logarithmic solution %%%%%%%%%%
uex = log(rad);