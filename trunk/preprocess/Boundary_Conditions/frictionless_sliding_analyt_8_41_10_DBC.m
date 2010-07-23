% DBCs

num_x = 41;
num_y = 10;

% left boundary
for i=1:num_y+1
  nodeID = i;
  analytic = fless_sliding_analyt_8_disp(x(nodeID),y(nodeID));
  dispbc(1,nodeID) = 1;
  ubar(1,nodeID) = analytic(1);
  dispbc(2,nodeID) = 1;
  ubar(2,nodeID) = analytic(2);
end;

% right boundary
for i=1:(num_y + 1)
  nodeID = num_x * (num_y + 1) + i;
  analytic = fless_sliding_analyt_8_disp(x(nodeID),y(nodeID));
  dispbc(1,nodeID) = 1;
  ubar(1,nodeID) = analytic(1);
  dispbc(2,nodeID) = 1;
  ubar(2,nodeID) = analytic(2);
end;

% dispbc(1,floor(num_x/2) * (num_y + 1) + ceil(num_y/2) + 1) = 1;
% dispbc(2,floor(num_x/2) * (num_y + 1) + ceil(num_y/2) + 1) = 1;
% 
% dispbc(1,ceil(num_x/2) * (num_y + 1) + ceil(num_y/2) + 1) = 1;
% dispbc(2,ceil(num_x/2) * (num_y + 1) + ceil(num_y/2) + 1) = 1;
% 
% dispbc(1,ceil(num_y/2) + 1) = 1;
% 
% dispbc(1,(num_x + 1) * (num_y + 1) - ceil(num_y/2)) = 1;
% 
% dispbc(1,(floor(num_x/2) - 1) * (num_y + 1) + ceil(num_y/2) + 1) = 1;
% dispbc(2,(floor(num_x/2) - 1) * (num_y + 1) + ceil(num_y/2) + 1) = 1;
% 
% dispbc(1,(ceil(num_x/2) + 1) * (num_y + 1) + ceil(num_y/2) + 1) = 1;
% dispbc(2,(ceil(num_x/2) + 1) * (num_y + 1) + ceil(num_y/2) + 1) = 1;