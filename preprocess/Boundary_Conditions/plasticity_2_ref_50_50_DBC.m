% DBCs

% mesh data
num_x = 20;
num_y = 20;

% left boundary
for i=1:(num_y + 1)
  dispbc(1,i) = 1;
  dispbc(2,i) = 1;
end;

