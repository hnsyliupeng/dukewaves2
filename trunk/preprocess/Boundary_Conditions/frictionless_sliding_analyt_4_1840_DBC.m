% DBCs

% bottom left corner
dispbc(1,6) = 1;
dispbc(2,6) = 1;

% bottom right corner
% dispbc(1,5) = 1;
dispbc(2,5) = 1;

% upper left corner
dispbc(1,3) = 1;
% dispbc(2,3) = 1;

% upper right corner
% dispbc(1,4) = 1;
% dispbc(2,4) = 1;

% upper cavity corner
dispbc(1,2) = 1;
% dispbc(2,2) = 1;

% lower cavity corner
dispbc(1,7) = 1;
% dispbc(2,7) = 1;

% upper left boundary
for i=8:23
  dispbc(1,i) = 1;
end;

% bottom left boundary
for i=98:113
  dispbc(1,i) = 1;
end;

% bottom boundary
for i=80:97
  dispbc(2,i) = 1;
end;