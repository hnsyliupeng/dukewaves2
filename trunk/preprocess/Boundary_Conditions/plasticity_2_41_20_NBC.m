% NBCs

% mesh data
num_x = 41;
num_y = 20;

% get nodes of right boundary
tractionnodes = [];
for i=1:(num_y + 1)
  nodeID = num_x * (num_y + 1) + i;
  tractionnodes = [tractionnodes nodeID];
end;

% define traction
FORCE(1) = struct('shape','','values',[],'nodes',[],'coords',[]);
FORCE(1).shape = 'constant';
FORCE(1).values = [1 0];
FORCE(1).nodes = tractionnodes;
