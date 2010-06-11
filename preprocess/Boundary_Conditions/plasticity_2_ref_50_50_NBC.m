% NBCs

% mesh data
num_x = 50;
num_y = 50;

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
