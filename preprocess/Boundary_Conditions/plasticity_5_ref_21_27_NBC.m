% NBC

% mesh data
num_x = 61;%21;%41;%61;%121;%161;
num_y = 61;%21;%41;%61;%121;%161;

p = -0.25e-3;         % load per length
loadlength = 4;   % length, that is loaded
load = p * loadlength / num_x; % nodal load

% bottom boundary
tractionnodes = [];
for i=1:(num_x + 1)
  nodeID = i * (num_y +1);
  tractionnodes = [tractionnodes nodeID];
end;

clear nodeID;

for nodeID = tractionnodes
  force(1,nodeID) = load;
end;
force(1,tractionnodes(1)) = force(1,tractionnodes(1)) / 2;
force(1,tractionnodes(end)) = force(1,tractionnodes(end)) / 2;

clear tractionnodes;

% % right boundary
% tractionnodes = [];
% for i=1:(num_y + 1)
%   nodeID = num_x * (num_y + 1) + i;
%   tractionnodes = [tractionnodes nodeID];
% end;
% 
% % define traction
% FORCE(1) = struct('shape','','values',[],'nodes',[],'coords',[]);
% FORCE(1).shape = 'constant';
% FORCE(1).nodes = tractionnodes;
% FORCE(1).values = [p/(length(tractionnodes) - 1) 0];
% 
% clear tractionnodes;


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
% FORCE(1).values = [p/num_x 0];
% 
% clear tractionnodes;



% % define traction
% FORCE(3) = struct('shape','','values',[],'nodes',[],'coords',[]);
% FORCE(3).shape = 'constant';
% FORCE(3).nodes = tractionnodes;
% FORCE(3).values = [0 0];
% 
% clear tractionnodes;
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