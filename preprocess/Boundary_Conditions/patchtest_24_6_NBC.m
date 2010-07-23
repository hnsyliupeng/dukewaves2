% NBCs

load = 10;%10;%0;
num_y = 6;

nodeload = load/num_y;

% forces in x-direction
force(1,169) = nodeload/2;
force(1,170) = nodeload;
force(1,171) = nodeload;
force(1,172) = nodeload;
force(1,173) = nodeload;
force(1,174) = nodeload;
force(1,175) = nodeload/2;



% % Input data for integrating the NBCs
% FORCE(1) = struct('shape','','values',[],'nodes',[],'coords',[]);
% FORCE(1).shape = 'constant';
% FORCE(1).values = [10 0];
% FORCE(1).nodes = 169:175;
