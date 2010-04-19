function [B] = sfderivatives(node,e,x,y,z)

xe = x(node(1:4,e));
ye = y(node(1:4,e));
ze = z(node(1:4,e));

%compute edge distances
x14 = xe(1) - xe(4);
x24 = xe(2) - xe(4);
x34 = xe(3) - xe(4);
y14 = ye(1) - ye(4);
y24 = ye(2) - ye(4);
y34 = ye(3) - ye(4);
z14 = ze(1) - ze(4);
z24 = ze(2) - ze(4);
z34 = ze(3) - ze(4);

%compute cofactor terms

cof11 = y24*z34 - y34*z24;
cof12 = y34*z14 - y14*z34;
cof13 = y14*z24 - y24*z14;
cof21 = z24*x34 - z34*x24;
cof22 = z34*x14 - z14*x34;
cof23 = z14*x24 - z24*x14;
cof31 = x24*y34 - x34*y24;
cof32 = x34*y14 - x14*y34;
cof33 = x14*y24 - x24*y14;

CC = [1,1,1,1; xe'; ye'; ze'];
Jcob = det(CC);

Cof = [cof11, cof12, cof13; cof21, cof22, cof23; cof31, cof32, cof33];

dN1dxdydz = (1/Jcob)*Cof*[1;0;0];
dN2dxdydz = (1/Jcob)*Cof*[0;1;0];
dN3dxdydz = (1/Jcob)*Cof*[0;0;1];
dN4dxdydz = (1/Jcob)*Cof*[-1;-1;-1];

%get sf derivatives
dNdx = -[dN1dxdydz(1) dN2dxdydz(1) dN3dxdydz(1) dN4dxdydz(1)];
dNdy = -[dN1dxdydz(2) dN2dxdydz(2) dN3dxdydz(2) dN4dxdydz(2)];
dNdz = -[dN1dxdydz(3) dN2dxdydz(3) dN3dxdydz(3) dN4dxdydz(3)];

B = [dNdx; dNdy; dNdz];