% DBCs

% fixed nodes
% fixed in x-direction
for i=1:19
    dispbc(1,i) = 1;
%     dispbc(2,i) = 1;
end;

% fixed in y-direction
dispbc(2,9) = 1;
dispbc(2,1378) = 1;

% right boundary
displacement = 0.04;
for i=1369:1387
    dispbc(1,i) = 1;
    ubar(1,i) = displacement;
end;