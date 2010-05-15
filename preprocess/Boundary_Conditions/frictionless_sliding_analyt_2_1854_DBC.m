% DBCs

% upper left corner
% dispbc(1,1) = 1;
% dispbc(2,1) = 1;

% upper right corner
% dispbc(1,2) = 1;
% dispbc(2,2) = 1;

% bottom right corner
% dispbc(1,3) = 1;
dispbc(2,3) = 1;

% bottom left corner
dispbc(1,4) = 1;
dispbc(2,4) = 1;

% left boundary
% for i=110:126
%   dispbc(1,i) = 1;
%   dispbc(2,i) = 1;
% end;

% bottom boundary
for i=66:109
%   dispbc(1,i) = 1;
  dispbc(2,i) = 1;
end;