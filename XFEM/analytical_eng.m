function [anal] = analytical_eng(x,y);

E = 1000.0;
pr = 0.3; 
P = -1;
h = 4;
I = h^3/12;
L = 16;
c = h/2;

anal(1) = -P*y*(L-x)/(E*I);

anal(2) = P*y*pr*(L-x)/(6*E*I);

anal(3) = 2*P*(1+pr)*(c^2 - y^2)/(2*E*I);

