% bc_conv9_NBC.m
%
% Boundary conditions for beam bending problem based on L = 16, set of
% structured mesh. This file is for Neumann BCs for 26 elements through the
% thickness and 104 elements over the length.
%

% Right hand side - parabolic distribution of downward shear, P = -1
force(2,2809) =   -51/35152;
force(2,2810) =  -149/17576;
force(2,2811) =  -287/17576;
force(2,2812) =  -413/17576;
force(2,2813) =  -527/17576;
force(2,2814) =  -629/17576;
force(2,2815) =  -719/17576;
force(2,2816) =  -797/17576;
force(2,2817) =  -863/17576;
force(2,2818) =  -917/17576;
force(2,2819) =  -959/17576;
force(2,2820) =  -989/17576;
force(2,2821) = -1007/17576;
force(2,2822) = -1013/17576;
force(2,2823) = -1007/17576;
force(2,2824) =  -989/17576;
force(2,2825) =  -959/17576;
force(2,2826) =  -917/17576;
force(2,2827) =  -863/17576;
force(2,2828) =  -797/17576;
force(2,2829) =  -719/17576;
force(2,2830) =  -629/17576;
force(2,2831) =  -527/17576;
force(2,2832) =  -413/17576;
force(2,2833) =  -287/17576;
force(2,2834) =  -149/17576;
force(2,2835) =   -51/35152;


% Left hand side - force boundary conditions

% y-forces to balance shear on RHS
force(2,1)  =   51/35152;
force(2,2)  =  149/17576;
force(2,3)  =  287/17576;
force(2,4)  =  413/17576;
force(2,5)  =  527/17576;
force(2,6)  =  629/17576;
force(2,7)  =  719/17576;
force(2,8)  =  797/17576;
force(2,9)  =  863/17576;
force(2,10) =  917/17576;
force(2,11) =  959/17576;
force(2,12) =  989/17576;
force(2,13) = 1007/17576;
force(2,15) = 1007/17576;
force(2,16) =  989/17576;
force(2,17) =  959/17576;
force(2,18) =  917/17576;
force(2,19) =  863/17576;
force(2,20) =  797/17576;
force(2,21) =  719/17576;
force(2,22) =  629/17576;
force(2,23) =  527/17576;
force(2,24) =  413/17576;
force(2,25) =  287/17576;
force(2,26) =  149/17576;
force(2,27) =   51/35152;

% x-forces to assure linear distribution of stresses
force(1,2)  = -144/169;
force(1,3)  = -132/169;
force(1,4)  = -120/169;
force(1,5)  = -108/169;
force(1,6)  =  -96/169;
force(1,7)  =  -84/169;
force(1,8)  =  -72/169;
force(1,9)  =  -60/169;
force(1,10) =  -48/169;
force(1,11) =  -36/169;
force(1,12) =  -24/169;
force(1,13) =  -12/169;
force(1,15) =   12/169;
force(1,16) =   24/169;
force(1,17) =   36/169;
force(1,18) =   48/169;
force(1,19) =   60/169;
force(1,20) =   72/169;
force(1,21) =   84/169;
force(1,22) =   96/169;
force(1,23) =  108/169;
force(1,24) =  120/169;
force(1,25) =  132/169;
force(1,26) =  144/169;