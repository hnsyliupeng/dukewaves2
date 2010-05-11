% NBCs

% define a linear traction in negative y-direction
FORCE(1) = struct('shape','','values',[],'nodes',[],'coords',[]);
FORCE(1).shape = 'linear';
FORCE(1).values = [0 0;0 -1];
FORCE(1).nodes = [1 81:132 4];

% FORCE(2) = struct('shape','','values',[],'nodes',[],'coords',[]);
% FORCE(2).shape = 'linear';
% FORCE(2).values = [0 0;0 1];
% FORCE(2).nodes = [2 68:-1:17 3];

% FORCE(1) = struct('shape','','values',[],'nodes',[],'coords',[]);
% FORCE(1).shape = 'constant';
% FORCE(1).values = [0 -1];
% FORCE(1).nodes = [1 81:132 4];


