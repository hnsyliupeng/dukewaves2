% DBCs

% x-direction
dispbc(1,1) = 1;
dispbc(1,2) = 1;
% dispbc(1,3) = 1;

% dispbc(1,22) = 1;
% dispbc(1,5) = 1;
% dispbc(1,4) = 1;


% y-direction
dispbc(2,2) = 1;
for i=7:20
  dispbc(2,i) = 1;
  dispbc(1,i) = 1;
end;

dispbc(2,3) = 1;