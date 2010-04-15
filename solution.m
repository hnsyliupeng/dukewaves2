clear
length = 40;
height = 4;

% number of elements in each direction
ndivl = 40;
ndivw = 16;

numele = ndivw*ndivl;
numnod = (ndivl+1)*(ndivw+1);


% set up nodal coordinates

for i = 1:(ndivl+1)
   for j=1:(ndivw+1)
      x((ndivw+1)*(i-1)+j) = (length/ndivl)*(i-1);
      y((ndivw+1)*(i-1)+j) = height/2 -(height/ndivw)*(j-1);
   end
end

% get analytical solution

for i = 1:numnod
    soln = analytical(x(i),y(i));
    disp(2*i-1)  = soln(1);
    disp(2*i)    = soln(2);
end


disp = disp*100;

for i = 1:numnod
x_def(i) = x(i) + disp(2*i-1);
y_def(i) = y(i) + disp(2*i);
end