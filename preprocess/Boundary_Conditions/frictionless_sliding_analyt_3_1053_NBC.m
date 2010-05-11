% NBCs

% define traction
FORCE(1) = struct('shape','','values',[],'nodes',[],'coords',[]);
FORCE(1).shape = 'constant';
FORCE(1).values = [0 -1];
FORCE(1).nodes = [5 202:266 2];