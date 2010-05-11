% DBCs

% dispbc(2,3) = 1;
dispbc(2,4) = 1;
dispbc(2,2) = 1;
dispbc(2,1) = 1;

dispbc(1,1) = 1;
dispbc(1,4) = 1;

for i=5:12
  dispbc(1,i) = 1;
  dispbc(2,i) = 1;
end;

for i=44:66
  dispbc(2,i) = 1;
end;

% for i=13:35
%   dispbc(2,i) = 1;
% end;