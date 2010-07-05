% NBCs

% mesh data
num_x = 125;
num_y = 5;

F = 44.48;%44.48;
h = 6.35;
fy = 3/4 * F / h / 12.7;

% get nodes of right boundary
tractionnodes = [];
for i=1:(num_y + 1)
  nodeID = num_x * (num_y + 1) + i;
  tractionnodes = [tractionnodes nodeID];
end;

% define traction
FORCE(1) = struct('shape','','values',[],'nodes',[],'coords',[]);
FORCE(1).shape = 'parabolic';
FORCE(1).nodes = tractionnodes;
FORCE(1).values = [0 fy/(length(tractionnodes) - 1)];

clear tractionnodes;

% % upper boundary
% tractionnodes = [];
% for i = 1:(num_x +1)
%   nodeID = i*(num_y + 1) - num_y;
%   tractionnodes = [tractionnodes nodeID];
% end;
% 
% % define traction
% FORCE(2) = struct('shape','','values',[],'nodes',[],'coords',[]);
% FORCE(2).shape = 'constant';
% FORCE(2).nodes = tractionnodes;
% FORCE(2).values = [0 -27.6/(length(tractionnodes)-1)];
% 
% clear tractionnodes;
% 
% % bottom boundary
% tractionnodes = [];
% for i = 1:(num_x +1)
%   nodeID = i*(num_y + 1);
%   tractionnodes = [tractionnodes nodeID];
% end;
% 
% % define traction
% FORCE(3) = struct('shape','','values',[],'nodes',[],'coords',[]);
% FORCE(3).shape = 'constant';
% FORCE(3).nodes = tractionnodes;
% FORCE(3).values = [0 27.6/(length(tractionnodes)-1)];

