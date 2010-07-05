% NBC

% mesh data
num_x = 147;
num_y = 49;

E = 1.0e+4;
p = 0;%1;%1.0e+4 * E;

% right boundary
tractionnodes = [];
for i=1:(num_y + 1)
  nodeID = num_x * (num_y + 1) + i;
  tractionnodes = [tractionnodes nodeID];
end;

% define traction
FORCE(1) = struct('shape','','values',[],'nodes',[],'coords',[]);
FORCE(1).shape = 'constant';
FORCE(1).nodes = tractionnodes;
FORCE(1).values = [p/(length(tractionnodes) - 1) 0];

clear tractionnodes;
% bottom boundary
tractionnodes = [];
for i=1:(num_x + 1)
  nodeID = i * (num_y +1);
  tractionnodes = [tractionnodes nodeID];
end;

% define traction
FORCE(2) = struct('shape','','values',[],'nodes',[],'coords',[]);
FORCE(2).shape = 'constant';
FORCE(2).nodes = tractionnodes;
FORCE(2).values = [0 -p*3/(length(tractionnodes)-1)];

clear tractionnodes;
