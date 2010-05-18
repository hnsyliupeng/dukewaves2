% NBCs

% define a linear traction
FORCE(1) = struct('shape','','values',[],'nodes',[],'coords',[]);
FORCE(1).shape = 'linear';
FORCE(1).values = [1 0;-1 0];
FORCE(1).nodes = [4:10 3];

FORCE(2) = struct('shape','','values',[],'nodes',[],'coords',[]);
FORCE(2).shape = 'linear';
FORCE(2).values = [-1 0;1 0];
FORCE(2).nodes = [2 41:46 1];

FORCE(3) = struct('shape','','values',[],'nodes',[],'coords',[]);
FORCE(3).shape = 'constant';
FORCE(3).values = [-2 0];
FORCE(3).nodes = [4:10 3];

FORCE(4) = struct('shape','','values',[],'nodes',[],'coords',[]);
FORCE(4).shape = 'constant';
FORCE(4).values = [2 0];
FORCE(4).nodes = [2 41:46 1];