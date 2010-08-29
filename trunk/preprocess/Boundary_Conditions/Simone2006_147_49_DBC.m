% DBC

% mesh data
num_x = 28;%28;%89;%147;
num_y = 19;%19;%57;%49;

% left boundary
for i=1:(num_y + 1)
  dispbc(1,i) = 1;
%   dispbc(2,i) = 1;
end;

% right boundary
for i=1:(num_y + 1)
  nodeID = num_x * (num_y + 1) + i;
  dispbc(1,nodeID) = 1;
  ubar(1,nodeID) = 3.0e-3;
end;

% bottom boundary
for i=1:(num_x + 1)
  nodeID = i * (num_y +1);
%   dispbc(1,nodeID) = 1;
  dispbc(2,nodeID) = 1;
%   ubar(2,nodeID) = -sqrt(3) * 1e-6;
end;


% upper boundary
for i = 1:(num_x +1)
  nodeID = i*(num_y + 1) - num_y;
%   dispbc(1,nodeID) = 1;
%   ubar(1,nodeID) = 0.03;%0.008;
  dispbc(2,nodeID) = 1;
%   ubar(2,nodeID) = -0.007;
end;