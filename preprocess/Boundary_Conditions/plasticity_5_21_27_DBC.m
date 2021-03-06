% DBC

% mesh data
num_x = 21;%21;%61;
num_y = 27;%27;%77;

% % left boundary
% for i=1:(num_y + 1)
%   dispbc(1,i) = 1;
% %   dispbc(2,i) = 1;
% end;
% 
% % right boundary
% for i=1:(num_y + 1)
%   nodeID = num_x * (num_y + 1) + i;
%   dispbc(1,nodeID) = 1;
% %   ubar(1,nodeID) = 0.0001;
% end;

% bottom boundary
for i=1:(num_x + 1)
  nodeID = i * (num_y +1);
  dispbc(1,nodeID) = 1;
  dispbc(2,nodeID) = 1;
%   ubar(2,nodeID) = -0.0001;
end;
% 
% dispbc(1,num_y + 1) = 1;


% % upper boundary
% for i = 1:(num_x +1)
%   nodeID = i*(num_y + 1) - num_y;
% %   dispbc(1,nodeID) = 1;
% %   ubar(1,nodeID) = 0.01;%0.008;
%   dispbc(2,nodeID) = 1;
% %   ubar(2,nodeID) = -0.1;
% end;

% move upper block
for i=1:(num_x+1)
  for j=1:(num_y - num_x - 1)
    nodeID = (i-1) * (num_y + 1) + j;
    dispbc(1,nodeID) = 1;
    ubar(1,nodeID) = 0.01;
    dispbc(2,nodeID) = 1;
  end;
end;