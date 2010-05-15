% Boundary Description

% left boundary
for i=1:IFnldivy
  boundary_nodes_left(1,i) = i;
  boundary_nodes_left(2,i) = i + 1;
end;

% right boundary
for i=1:IFnldivy
  boundary_nodes_right(1,i) = IFnldivx * (IFnldivy+1) + i;
  boundary_nodes_right(2,i) = IFnldivx * (IFnldivy+1) + i +1;
end;

% bottom boundary
for i=2:IFnldivx
  boundary_nodes_bottom(1,i-1) = (IFnldivy + 1) * i;
  boundary_nodes_bottom(2,i-1) = (IFnldivy + 1) * (i + 1);
end;

% top boundary
for i=1:IFnldivx-1
  boundary_nodes_top(1,i) = (IFnldivy + 1) * i + 1;
  boundary_nodes_top(2,i) = (IFnldivy + 1) * (i + 1) + 1;
end;

boundary_nodes = [boundary_nodes_left boundary_nodes_bottom ...
  boundary_nodes_right boundary_nodes_top];

numboundele = size(boundary_nodes,2);