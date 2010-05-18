% NBCs

num_x = 161;
num_y = 40;

% define a linear traction
% left boundary
FORCE(1) = struct('shape','','values',[],'nodes',[],'coords',[]);
FORCE(1).shape = 'linear';
FORCE(1).values = [1 0;-1 0];
FORCE(1).nodes = [1:num_y+1];
  
% right boundary
FORCE(2) = struct('shape','','values',[],'nodes',[],'coords',[]);
FORCE(2).shape = 'linear';
FORCE(2).values = [-1 0;1 0];
FORCE(2).nodes = [num_x * (num_y + 1) + 1:(num_x + 1) * (num_y + 1)];

% load = 1/10;
% 
% % bending moments
% for i=1:7
%   % right end
%   force(1,1142-i) = load * i;
%   force(1,1141+i) = -load * i;
%   
%   % left end
%   force(1,8-i) = -load * i;
%   force(1,7+i) = load * i;
% end;