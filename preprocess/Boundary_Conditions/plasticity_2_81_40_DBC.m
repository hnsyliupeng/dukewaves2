% DBCs

% mesh data
num_x = 81;
num_y = 40;

% left boundary
for i=1:(num_y + 1)
  dispbc(1,i) = 1;
  dispbc(2,i) = 1;
end;

% right boundary
displacement = 0.01;%0.5-0.022191764346110;%0.5;%0.001;
for i=1:(num_y + 1)
  nodeID = num_x * (num_y + 1) + i;
  dispbc(1,nodeID) = 1;
  dispbc(2,nodeID) = 1;
  ubar(2,nodeID) = displacement;
end;

% nodeID = (num_x+1)*(num_y+1);
% dispbc(1,nodeID) = 1;
% dispbc(2,nodeID) = 1;
% 
% nodeID = nodeID - 1;
% dispbc(1,nodeID) = 1;
% dispbc(2,nodeID) = 1;
% 
% nodeID = nodeID - 1;
% dispbc(1,nodeID) = 1;
% dispbc(2,nodeID) = 1;



% nodeID = (num_x)*(num_y+1)+1;
% dispbc(1,nodeID) = 1;
% dispbc(2,nodeID) = 1;