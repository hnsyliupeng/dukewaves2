% NBC

% mesh data
num_x = 147;
num_y = 49;

E = 1.0e+4;
p = 1;%1.0e+4 * E;

% right boundary
tractionnodes = [];
for i=1:(num_y + 1)
  nodeID = num_x * (num_y + 1) + i;
  tractionnodes = [tractionnodes nodeID];
end;

% define traction
FORCE(1) = struct('shape','','values',[],'nodes',[],'coords',[]);
FORCE(1).shape = 'constant';
FORCE(1).nodes = tractionnodes;
FORCE(1).values = [p/(length(tractionnodes) - 1) 0];

clear tractionnodes;


% % bottom boundary
% tractionnodes = [];
% for i=1:(num_x + 1)
%   nodeID = i * (num_y +1);
%   tractionnodes = [tractionnodes nodeID];
% end;
% 
% % define traction
% FORCE(1) = struct('shape','','values',[],'nodes',[],'coords',[]);
% FORCE(1).shape = 'constant';
% FORCE(1).nodes = tractionnodes;
% FORCE(1).values = [-44.1*20/(length(tractionnodes)-1) 0];
% 
% clear tractionnodes;
% 
% 
% % % upper boundary
% % tractionnodes = [];
% % for i = 1:(num_x +1)
% %   nodeID = i*(num_y + 1) - num_y;
% %   tractionnodes = [tractionnodes nodeID];
% % end;
% % 
% % % define traction
% % FORCE(3) = struct('shape','','values',[],'nodes',[],'coords',[]);
% % FORCE(3).shape = 'constant';
% % FORCE(3).nodes = tractionnodes;
% % FORCE(3).values = [0 0];
% % 
% % clear tractionnodes;
% % 
% % 
% % % left boundary
% % tractionnodes = [];
% % for i = 1:(num_x +1)
% %   nodeID = i;
% %   tractionnodes = [tractionnodes nodeID];
% % end;
% % 
% % % define traction
% % FORCE(4) = struct('shape','','values',[],'nodes',[],'coords',[]);
% % FORCE(4).shape = 'constant';
% % FORCE(4).nodes = tractionnodes;
% % FORCE(4).values = [1*4/(length(tractionnodes) - 1) 0];
% % 
% % clear tractionnodes;