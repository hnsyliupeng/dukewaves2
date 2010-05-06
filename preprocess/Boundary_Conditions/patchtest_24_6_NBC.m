% NBCs

load=10/6;

% forces in x-direction
force(1,169) = load/2;
force(1,170) = load;
force(1,171) = load;
force(1,172) = load;
force(1,173) = load;
force(1,174) = load;
force(1,175) = load/2;

% forces in y-direction


% % Input data for integrating the NBCs
% FORCE(1) = struct('shape','','values',[],'nodes',[],'coords',[]);
% FORCE(1).shape = 'constant';
% FORCE(1).values = [10 0];
% FORCE(1).nodes = 169:175;
