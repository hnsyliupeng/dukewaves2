% NBCs

% mesh data
num_x = 20;
num_y = 20;

% get nodes of right boundary
tractionnodes = [];
for i=1:(num_y + 1)
  nodeID = num_x * (num_y + 1) + i;
  tractionnodes = [tractionnodes nodeID];
end;

fy = 0.03 * 4 / (length(tractionnodes) - 1);

% define traction
FORCE(1) = struct('shape','','values',[],'nodes',[],'coords',[]);
FORCE(1).shape = 'constant';
FORCE(1).values = [0 fy];
FORCE(1).nodes = tractionnodes;

% define traction
FORCE(2) = struct('shape','','values',[],'nodes',[],'coords',[]);
FORCE(2).shape = 'linear';
FORCE(2).values = [1.5 0; -1.5 0];
FORCE(2).nodes = tractionnodes;

% FORCE(2).values = FORCE(2).values ./ (length(tractionnodes)-1);
