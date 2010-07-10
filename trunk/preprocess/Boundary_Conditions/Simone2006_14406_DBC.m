% DBC

% left boundary
leftboundarynodes = [4 267:324 1];
for i=leftboundarynodes
  nodeID = i;
  dispbc(1,nodeID) = 1;
%   dispbc(2,i) = 1;
end;

% right boundary
rightboundarynodes = [2 107:164 3];
for i=rightboundarynodes
  nodeID = i;
  dispbc(1,nodeID) = 1;
  ubar(1,nodeID) = 3.0e-6;
end;

% bottom boundary
bottomboundarynodes = [3 165:266 4];
for i=bottomboundarynodes
  nodeID = i;
  %   dispbc(1,nodeID) = 1;
  dispbc(2,nodeID) = 1;
%   ubar(2,nodeID) = -0.0001;
end;


% upper boundary
upperboundarynodes = [1 5:106 2];
for i=upperboundarynodes
  nodeID = i;
%   dispbc(1,nodeID) = 1;
%   ubar(1,nodeID) = 0.03;%0.008;
  dispbc(2,nodeID) = 1;
%   ubar(2,nodeID) = -0.007;
end;