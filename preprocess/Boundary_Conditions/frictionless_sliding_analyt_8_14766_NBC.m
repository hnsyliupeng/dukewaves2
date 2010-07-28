% NBCs

% define a linear traction
% right boundary
FORCE(1) = struct('shape','','values',[],'nodes',[],'coords',[]);
FORCE(1).shape = 'linear';
FORCE(1).values = [-1 0;1 0];
FORCE(1).nodes = [5 241:258 6 276:-1:259 7];

% left boundary
FORCE(2) = struct('shape','','values',[],'nodes',[],'coords',[]);
FORCE(2).shape = 'linear';
FORCE(2).values = [1 0;-1 0];
FORCE(2).nodes = [10 13:30 1 31:48 2];

% % right boundary
% FORCE(3) = struct('shape','','values',[],'nodes',[],'coords',[]);
% FORCE(3).shape = 'constant';
% FORCE(3).values = [-1 0];
% FORCE(3).nodes = [4:22 3];
% 
% % left boundary
% FORCE(4) = struct('shape','','values',[],'nodes',[],'coords',[]);
% FORCE(4).shape = 'constant';
% FORCE(4).values = [1 0];
% FORCE(4).nodes = [2 101:118 1];