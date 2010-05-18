% DBCs

num_x = 321;
num_y = 80;

dispbc(1,floor(num_x/2) * (num_y + 1) + ceil(num_y/2) + 1) = 1;
dispbc(2,floor(num_x/2) * (num_y + 1) + ceil(num_y/2) + 1) = 1;

dispbc(1,ceil(num_x/2) * (num_y + 1) + ceil(num_y/2) + 1) = 1;
dispbc(2,ceil(num_x/2) * (num_y + 1) + ceil(num_y/2) + 1) = 1;

dispbc(1,ceil(num_y/2) + 1) = 1;

dispbc(1,(num_x + 1) * (num_y + 1) - ceil(num_y/2)) = 1;

dispbc(1,(floor(num_x/2) - 1) * (num_y + 1) + ceil(num_y/2) + 1) = 1;
dispbc(2,(floor(num_x/2) - 1) * (num_y + 1) + ceil(num_y/2) + 1) = 1;

dispbc(1,(ceil(num_x/2) + 1) * (num_y + 1) + ceil(num_y/2) + 1) = 1;
dispbc(2,(ceil(num_x/2) + 1) * (num_y + 1) + ceil(num_y/2) + 1) = 1;