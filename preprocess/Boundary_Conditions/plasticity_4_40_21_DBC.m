% DBC

% mesh data
num_x = 40;
num_y = 21;

% % bottom boundary
% for i=1:(num_x + 1)
%   nodeID = i * (num_y +1);
%   dispbc(1,nodeID) = 1;
%   dispbc(2,nodeID) = 1;
% end;

% upper boundary
for i = 1:(num_x + 1)
  nodeID = i*(num_y + 1) - num_y;
  dispbc2(1,nodeID) = 1;
  ubar2(1,nodeID) = 0.008;
  dispbc(2,nodeID) = 1;
  ubar(2,nodeID) = -0.007;
  dispbc3(1,nodeID) = 1;
%   ubar3(2,nodeID) = -0.007;
end;

% fix the bottom block
for i=1:10
  for j=0:(num_x)
    nodeID = j*22 + i + 12;
    dispbc(1,nodeID) = 1;
    dispbc(2,nodeID) = 1;
  end;
end;