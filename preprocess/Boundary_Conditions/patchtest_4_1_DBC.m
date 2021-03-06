% DBCs

% fixed nodes
% fixed in x-direction
dispbc(1,1) = 1;
dispbc(1,2) = 1;

% fixed in y-direction
dispbc(2,1) = 1;
dispbc(2,2) = 1;

% right boundary
for i=9:10
  dispbc(2,i) = 1;
  ubar(2,i) = -0.5;
end;