% DBC

% mesh data
num_x = 80;
num_y = 41;

% % upper boundary
% for i = 1:(num_x +1)
%   nodeID = i*(num_y + 1) - num_y;
%   dispbc(1,nodeID) = 1;
%   ubar(1,nodeID) = 0.03;%0.008;
%   dispbc(2,nodeID) = 1;
%   ubar(2,nodeID) = -0.007;
%   dispbc2(1,nodeID) = 1;
%   ubar2(1,nodeID) = 0.008;
% end;

% fix the bottom block
for i=1:20
  for j=0:(num_x)
    nodeID = j*42 + i + 22;
    dispbc(1,nodeID) = 1;
    dispbc(2,nodeID) = 1;
  end;
end;