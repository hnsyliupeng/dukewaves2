% NBCs

for i=0:81
  forcenodes(i+1) = i * (6 + 1) + 1; 
end;

% define a linear traction
FORCE(1) = struct('shape','','values',[],'nodes',[],'coords',[]);
FORCE(1).shape = 'constant';
FORCE(1).values = [0 -0.1];
FORCE(1).nodes = forcenodes;




% load = -1/82;
% 
% for i=0:81
%   force(2,7*i + 1) = load/10; 
% end;