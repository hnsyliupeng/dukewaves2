% NBCs

load=10/98;

% forces in x-direction
force(1,38809) = load/2;

for i=38810:38906
    force(1,i) = load;
end;

force(1,38907) = load/2;

% forces in y-direction
