% DBCs

num_x = 750;
num_y = 31;

% left boundary
for i=1:(num_y + 1)
  dispbc(1,i) = 1;
%   dispbc(2,i) = 1;
end;


dispbc(2,1) = 1;
dispbc(2,num_y+1) = 1;

dispbc(2,3) = 1;
dispbc(2,4) = 1;
