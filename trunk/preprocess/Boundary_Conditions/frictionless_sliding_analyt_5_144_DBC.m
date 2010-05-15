% DBCs

% bottom left corner
dispbc(1,2) = 1;
dispbc(2,2) = 1;

% bottom right corner
dispbc(1,3) = 1;
dispbc(2,3) = 1;

% upper left corner
dispbc(1,1) = 1;
dispbc(2,1) = 1;

% upper right corner
dispbc(1,4) = 1;
dispbc(2,4) = 1;

% left boundary
for i=5:12
  dispbc(1,i) = 1;
  dispbc(2,i) = 1;
end;

% right boundary
for i=71:78
  dispbc(1,i) = 1;
  dispbc(2,i) = 1;
end;

% for i=79:87
%   dispbc(2,i) = 1;
% end;
% for i=62:70
%   dispbc(2,i) = 1;
% end;
