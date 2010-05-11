% DBCs

% bottom left corner
dispbc(1,4) = 1;
dispbc(2,4) = 1;

% bottom right corner
dispbc(1,3) = 1;
dispbc(2,3) = 1;

% upper left corner
dispbc(1,5) = 1;
% dispbc(2,5) = 1;

% upper right corner
dispbc(1,2) = 1;
% dispbc(2,2) = 1;

% bottom boundary
for i=72:136
  dispbc(2,i) = 1;
end;