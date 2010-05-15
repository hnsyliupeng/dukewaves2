% NBCs

% define a linear traction
FORCE(1) = struct('shape','','values',[],'nodes',[],'coords',[]);
FORCE(1).shape = 'linear';
FORCE(1).values = [0 0;0 -1];
FORCE(1).nodes = [1 5:74 2];