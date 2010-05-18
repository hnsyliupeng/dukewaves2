% NBCs

num_x = 641;
num_y = 48;

for i=0:num_x
  forcenodes(i+1) = i * (num_y + 1) + 1; 
end;

% define a linear traction
FORCE(1) = struct('shape','','values',[],'nodes',[],'coords',[]);
FORCE(1).shape = 'constant';
FORCE(1).values = [0 -0.1];
FORCE(1).nodes = forcenodes;

