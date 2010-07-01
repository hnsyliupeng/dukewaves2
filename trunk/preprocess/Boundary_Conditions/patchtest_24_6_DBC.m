% DBCs

% fixed nodes
% fixed in x-direction
dispbc(1,1) = 1;
dispbc(1,2) = 1;
dispbc(1,3) = 1;
dispbc(1,4) = 1;
dispbc(1,5) = 1;
dispbc(1,6) = 1;
dispbc(1,7) = 1;

% fixed in y-direction
% dispbc(2,1) = 1;
% dispbc(2,2) = 1;
% dispbc(2,3) = 1;
dispbc(2,4) = 1;
% dispbc(2,5) = 1;
% dispbc(2,6) = 1;
% dispbc(2,7) = 1;

% right boundary with prescribed displacement
displacement = 0.04;
for i=169:175
  dispbc(1,i) = 1;
  ubar(1,i) = displacement;
end;
