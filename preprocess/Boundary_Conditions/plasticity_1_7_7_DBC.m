% DBCs

% mesh
num_x = 7;
num_y = 7;

% bottom boundary
for i=1:(num_x + 1)
  nodeID = (num_y + 1) * i;
  dispbc(1,nodeID) = 1;
  dispbc(2,nodeID) = 1;
end;

% left boundary
for i=1:(num_y + 1)
  nodeID = i;
  dispbc(1,nodeID) = 1;
  dispbc(2,nodeID) = 1;
end;

% % release the DBCs on enriched nodes
% dispbc(1,4) = 0;
% dispbc(2,4) = 0;
% dispbc(1,5) = 0;
% dispbc(2,5) = 0;