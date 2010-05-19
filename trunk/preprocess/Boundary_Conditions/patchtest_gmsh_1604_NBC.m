% NBCs

% Input data for integrating the NBCs
FORCE(1) = struct('shape','','values',[],'nodes',[],'coords',[]);
FORCE(1).shape = 'constant';
FORCE(1).nodes = [4:16 3];
FORCE(1).values = [-1/13 0];


% FORCE(1) = struct('shape','','values',[],'nodes',[],'coords',[]);
% FORCE(1).shape = 'linear';
% FORCE(1).values = [10 0; 0 0];
% FORCE(1).nodes = [4 5 6 7 8 9 10 11 12 13 14 15 16 3];