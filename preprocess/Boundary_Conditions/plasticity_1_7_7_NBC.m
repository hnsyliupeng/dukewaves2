% NBCs

% mesh data
num_x = 7;
num_y = 7;

% get nodes of upper boundary
tractionnodes = [];
for i=1:(num_x + 1)
  nodeID = i * (num_y + 1) - num_y;
  tractionnodes = [tractionnodes nodeID];
end;

% define traction on upper boundary
FORCE(1) = struct('shape','','values',[],'nodes',[],'coords',[]);
FORCE(1).shape = 'constant';
FORCE(1).values = [-10 0];
FORCE(1).nodes = tractionnodes;

clear tractionnodes;

% get nodes of right boundary
tractionnodes = [];
for i=1:(num_y + 1)
  nodeID = num_x * (num_y + 1) + i;
  tractionnodes = [tractionnodes nodeID];
end;

% define traction on right boundary
FORCE(2) = struct('shape','','values',[],'nodes',[],'coords',[]);
FORCE(2).shape = 'constant';
FORCE(2).values = [0 -10];
FORCE(2).nodes = tractionnodes;