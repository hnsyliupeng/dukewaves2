% NBCs

FORCE(1) = struct('shape','','values',[],'nodes',[],'coords',[]);
FORCE(1).shape = 'constant';
FORCE(1).nodes = [4:17 3];
FORCE(1).values = [-0.25*4/(length(FORCE(1).nodes) - 1)*10 0];


