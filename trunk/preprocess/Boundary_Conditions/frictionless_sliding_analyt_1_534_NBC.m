% NBCs

% define a linear traction in negative y-direction
FORCE(1) = struct('shape','','values',[],'nodes',[],'coords',[]);
FORCE(1).shape = 'linear';
FORCE(1).values = [0 0;0 -1];
FORCE(1).nodes = [1 47:76 4];
% 
% FORCE(2) = struct('shape','','values',[],'nodes',[],'coords',[]);
% FORCE(2).shape = 'linear';
% FORCE(2).values = [1 0;0 0];
% FORCE(2).nodes = [4:10 3];

% FORCE(1) = struct('shape','','values',[],'nodes',[],'coords',[]);
% FORCE(1).shape = 'constant';
% FORCE(1).values = [0 -10];
% FORCE(1).nodes = [1 47:76 4];
% 
% FORCE(2) = struct('shape','','values',[],'nodes',[],'coords',[]);
% FORCE(2).shape = 'constant';
% FORCE(2).values = [-10/4 0];
% FORCE(2).nodes = [4:10 3];



