% NBCs

load = -1;

nodenum = 27;

nodevec = [3 29:53 2];

for i = 1:length(nodevec)
  phi = pi/2 * (i-1) / (length(nodevec)-1);
  force(1,nodevec(i)) = load * cos(phi);
  force(2,nodevec(i)) = load * sin(phi);
end;

force(2,3) = force(2,3) / 2;
force(1,2) = force(1,2) / 2;