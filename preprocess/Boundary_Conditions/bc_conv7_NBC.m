% bc_conv7_NBC.m
%
% Boundary conditions for beam bending problem based on L = 16, set of
% structured mesh.  This file is for Neumann BCs 18 elements through the 
% thickness and 72 elements over the length.
%

% Right hand side - parabolic distribution of downward shear, P = -1
force(2,1369) =  -35/11664;
force(2,1370) = -101/5832;
force(2,1371) = -191/5832;
force(2,1372) = -269/5832;
force(2,1373) = -335/5832;
force(2,1374) = -389/5832;
force(2,1375) = -431/5832;
force(2,1376) = -461/5832;
force(2,1377) = -479/5832;
force(2,1378) = -485/5832;
force(2,1379) = -479/5832;
force(2,1380) = -461/5832;
force(2,1381) = -431/5832;
force(2,1382) = -389/5832;
force(2,1383) = -335/5832;
force(2,1384) = -269/5832;
force(2,1385) = -191/5832;
force(2,1386) = -101/5832;
force(2,1387) =  -35/11664;


% Left hand side - force boundary conditions

% y-forces to balance shear on RHS
force(2,1)  =  35/11664;
force(2,2)  = 101/5832;
force(2,3)  = 191/5832;
force(2,4)  = 269/5832;
force(2,5)  = 335/5832;
force(2,6)  = 389/5832;
force(2,7)  = 431/5832;
force(2,8)  = 461/5832;
force(2,9)  = 479/5832;
force(2,11) = 479/5832;
force(2,12) = 461/5832;
force(2,13) = 431/5832;
force(2,14) = 389/5832;
force(2,15) = 335/5832;
force(2,16) = 269/5832;
force(2,17) = 191/5832;
force(2,18) = 101/5832;
force(2,19) =  35/11664;

% x-forces to assure linear distribution of stresses
force(1,2)  = -32/27;
force(1,3)  = -28/27;
force(1,4)  = -24/27;
force(1,5)  = -20/27;
force(1,6)  = -16/27;
force(1,7)  = -12/27;
force(1,8)  =  -8/27;
force(1,9)  =  -4/27;
force(1,11) =   4/27;
force(1,12) =   8/27;
force(1,13) =  12/27;
force(1,14) =  16/27;
force(1,15) =  20/27;
force(1,16) =  24/27;
force(1,17) =  28/27;
force(1,18) =  32/27;