% DBC

% mesh data
num_x = 40;
num_y = 21;

% % left boundary
% for i=1:(num_y + 1)
%   dispbc(1,i) = 1;
%   dispbc(2,i) = 1;
% end;

% dispbc(2,5) = 1;

% dispbc(1,1) = 1;
% dispbc(1,num_y + 1) = 1;
% dispbc(2,num_y + 1) = 1;

% % right boundary
% for i=1:(num_y + 1)
%   nodeID = num_x * (num_y + 1) + i;
%   dispbc(2,nodeID) = 1;
%   ubar(2,nodeID) = 0.0001;
% end;

% bottom boundary
for i=1:(num_x + 1)
  nodeID = i * (num_y +1);
  dispbc(1,nodeID) = 1;
  dispbc(2,nodeID) = 1;
end;
% dispbc(1,55) = 1;

% % upper boundary
% for i = 1:(num_x +1)
%   nodeID = i*(num_y + 1) - num_y;
% %   dispbc(1,nodeID) = 1;
% %   ubar(1,nodeID) = 0.03;%0.008;
%   dispbc(2,nodeID) = 1;
%   ubar(2,nodeID) = -0.007;
% end;
% dispbc(1,881) = 1;
% ubar(1,881) = 0.008;

% fix the bottom block
for i=1:10
  for j=0:(num_x)
    nodeID = j*22 + i + 12;
    dispbc(1,nodeID) = 1;
    dispbc(2,nodeID) = 1;
  end;
end;