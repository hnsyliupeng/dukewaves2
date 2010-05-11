% DBCs

% x-direction
dispbc(1,1) = 1;
dispbc(1,2) = 1;
% dispbc(1,3) = 1;

% dispbc(1,73) = 1;
% dispbc(1,76) = 1;
% dispbc(2,73) = 1;
% dispbc(2,76) = 1;

% for i=3:16
%   dispbc(1,i) = 1;
% end;
for i=69:80
  dispbc(1,i) = 1;
end;
% for i=17:68
%   dispbc(2,i) = 1;
% end;
for i=81:132
  dispbc(2,i) = 1;
end;

dispbc(2,4) = 1;

dispbc(2,1) = 1;
% dispbc(2,3) = 1;
% dispbc(2,2) = 1;
% y-direction
