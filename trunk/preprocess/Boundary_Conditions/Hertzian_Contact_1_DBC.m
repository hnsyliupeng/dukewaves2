% DBCs

% enriched nodes with DBCs: 57 and 58

% fixed nodes
% fixed in x-direction
for i=32:66
    dispbc(1,i) = 1;
    ubar(1,i) = 0;
end;

dispbc(1,1) = 1;
ubar(1,1) = 0;

dispbc(1,5) = 1;
ubar(1,5) = 0;

% fixed in y-direction
dispbc(2,1) = 1;
ubar(2,1) = 0;

for i=67:81
    dispbc(2,i) = 1;
    ubar(2,i) = 0;
end;

dispbc(2,2) = 1;
ubar(2,2) = 0;




