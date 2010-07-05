% NBC

% boundary nodes
leftboundarynodes = [4 267:324 1];
rightboundarynodes = [2 107:164 3];
bottomboundarynodes = [3 165:266 4];
upperboundarynodes = [1 5:106 2];

E = 1.0e+4;
p = 0;%1;%1.0e+4 * E;

% right boundary
% define traction
FORCE(1) = struct('shape','','values',[],'nodes',[],'coords',[]);
FORCE(1).shape = 'constant';
FORCE(1).nodes = rightboundarynodes;
FORCE(1).values = [p/(length(tractionnodes) - 1) 0];

clear tractionnodes;


% bottom boundary
% define traction
FORCE(2) = struct('shape','','values',[],'nodes',[],'coords',[]);
FORCE(2).shape = 'constant';
FORCE(2).nodes = bottomboundarynodes;
FORCE(2).values = [0 -p*3/(length(tractionnodes)-1)];

clear tractionnodes;
