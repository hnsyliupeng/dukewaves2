function [volumep] = l2volumelszero(node,e,x,y,z)

xe = x(node(1:4,e));
ye = y(node(1:4,e));
ze = z(node(1:4,e));
nlink = 4;
volumep = tetraunitvolume('DUMMY')*parallelipipedvolume3d(xe, ye, ze);
