% DBC

% left upper corner
dispbc(1,3) = 1;
% dispbc(2,3) = 1;

% left lower corner
dispbc(1,4) = 1;
% dispbc(2,4) = 1;

% bottom left corner
% dispbc(1,5) = 1;
dispbc(2,5) = 1;

% bottom right corner
% dispbc(1,2) = 1;
dispbc(2,2) = 1;

% left boundary
for i=6:7
  dispbc(1,i) = 1;
%   dispbc(2,i) = 1;
end;

% right boundary
for i=8:9
%   dispbc(1,i) = 1;
  dispbc(2,i) = 1;
end;