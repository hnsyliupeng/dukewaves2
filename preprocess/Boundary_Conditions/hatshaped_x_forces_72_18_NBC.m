% NBCs

load=10/18;

% forces in x-direction
for i=0:9
    force(1,1369+i) = 3*i*load;
    force(1,1387-i) = 3*i*load;
end;

% forces in y-direction
