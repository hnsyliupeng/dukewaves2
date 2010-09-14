% DBC

% dy = 0.001;
% dx = -0.3 * dy;
dx = 0.001;
dy = dx;

% mesh data
num_x = 50;%40;
num_y = 50;%40;

% left boundary
for i=1:(num_y + 1)
  nodeID = i;
%   dispbc(1,nodeID) = 1;
%   ubar(1,nodeID) = -dx/2;
  dispbc(2,nodeID) = 1;
  ubar(2,nodeID) = -dy/2;
end;

% right boundary
for i=1:(num_y + 1)
  nodeID = num_x * (num_y + 1) + i;
%   dispbc(1,nodeID) = 1;
%   ubar(1,nodeID) = dx/2;
  dispbc(2,nodeID) = 1;
  ubar(2,nodeID) = dy/2;
end;

% bottom boundary
for i=1:(num_x + 1)
  nodeID = i * (num_y +1);
  dispbc(1,nodeID) = 1;
  ubar(1,nodeID) = -dx/2;
%   dispbc(2,nodeID) = 1;
%   ubar(2,nodeID) = -dy/2;
end;


% upper boundary
for i = 1:(num_x +1)
  nodeID = i*(num_y + 1) - num_y;
  dispbc(1,nodeID) = 1;
  ubar(1,nodeID) = dx/2;
%   dispbc(2,nodeID) = 1;
%   ubar(2,nodeID) = dy/2;
end;