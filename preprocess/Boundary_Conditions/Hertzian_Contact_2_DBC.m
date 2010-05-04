% DBCs

% fixed nodes
% fixed in x-direction
for i=43:90
    dispbc(1,i) = 1;
    ubar(1,i) = 0;
end;

dispbc(1,1) = 1;
ubar(1,1) = 0;

dispbc(1,2) = 1;
ubar(1,2) = 0;

dispbc(1,160) = 1;
ubar(1,160) = 0;

dispbc(1,161) = 1;
ubar(1,16) = 0;

% fixed in y-direction
dispbc(2,1) = 1;
ubar(2,1) = 0;

for i=43:90
    dispbc(2,i) = 1;
    ubar(2,i) = 0;
end;

dispbc(2,2) = 1;
ubar(2,2) = 0;




0