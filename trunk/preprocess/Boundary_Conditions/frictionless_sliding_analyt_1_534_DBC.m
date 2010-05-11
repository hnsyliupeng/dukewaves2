% DBCs

% x-direction
% dispbc(1,1) = 1;
dispbc(1,2) = 1;
% dispbc(1,3) = 1;

% for i=41:46
%   dispbc(1,i) = 1;
% end;
% for i=3:10
%   dispbc(1,i) = 1;
% end;


% y-direction
dispbc(2,2) = 1;
for i=11:40
  dispbc(2,i) = 1;
  dispbc(1,i) = 1;
end;

dispbc(2,3) = 1;