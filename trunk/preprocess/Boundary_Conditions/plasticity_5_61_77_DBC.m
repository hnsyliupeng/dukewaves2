% DBC

% mesh data
num_x = 141;%11;%21;%41;%61;%81;%141;
num_y = 176;%13;%27;%51;%77;%101;%176;

% upper boundary
for i = 1:(num_x +1)
  nodeID = i*(num_y + 1) - num_y;
  dispbc(1,nodeID) = 1;
  ubar(1,nodeID) = 0.005;
  dispbc(2,nodeID) = 1;  
end;

% fix bottom block
for i=1:(num_x+1)
  for j=1:(num_y - num_x)
    nodeID = (num_x + 1) + j + (i-1) * (num_y + 1);
    dispbc(1,nodeID) = 1;
    dispbc(2,nodeID) = 1;
  end;
end;