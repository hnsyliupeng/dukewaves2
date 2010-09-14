% frictionless_sliding1_8_2_DBC.m

% left boundary
% x-displacement
dispbc(1,1) = 1;
dispbc(1,2) = 1;
dispbc(1,3) = 1;

% y-displacement
dispbc(2,2) = 1;

% right boundary
% x- displacement
dx = -0.004;
dispbc(1,25) = 1;
ubar(1,25) = dx;
dispbc(1,26) = 1;
ubar(1,26) = dx;
dispbc(1,27) = 1;
ubar(1,27) = dx;

% y-displacement
dispbc(2,26) = 1;