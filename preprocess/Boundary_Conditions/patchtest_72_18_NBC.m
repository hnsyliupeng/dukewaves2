% NBCs

load=0;%10/18;

% forces in x-direction
force(1,1369) = load/2;

for i=1370:1386
    force(1,i) = load;
end;

force(1,1387) = load/2;

% forces in y-direction
