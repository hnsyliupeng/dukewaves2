% NBCs

load = 1/10;

% bending moments
for i=1:7
  % right end
  force(1,1142-i) = load * i;
  force(1,1141+i) = -load * i;
  
  % left end
  force(1,8-i) = -load * i;
  force(1,7+i) = load * i;
end;