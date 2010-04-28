% NBCs

load = 10;

% Right hand side - parabolic distribution

for i=1:9
    force(1,1378+i) = (1-(i^2)/81)*load;
    force(1,1378-i) = (1-(i^2)/81)*load;
end;
force(1,1378) = load;