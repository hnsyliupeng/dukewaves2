% NBCs

load = 10;

% Right hand side
% force(1,2) = load/2;
% force(1,3) = load/2;
% force(1,5) = load;
% force(1,6) = load;
% force(1,7) = load;

% define a structure, that defines the properties of the applied forces
% first set of tractions
% nodes_with_force = [3 5 6 7 2];
% forceshape = 'parabolic'; 
% values = [0 0; 10 0];        % [x-force y-force]
% FORCE(1) = struct('shape',forceshape,'values',values, ...
%     'nodes',nodes_with_force,'coords',[]);

% nodes_with_force = [4 14 15 16 3];
% forceshape = 'linear'; 
% values = [0 0; 0 -100];        % [x-force y-force]
% FORCE(2) = struct('shape',forceshape,'values',values, ...
%     'nodes',nodes_with_force,'coords',[]);
% 
% nodes_with_force = [1 10 9 8 2];
% forceshape = 'linear'; 
% values = [0 0; 0 100];        % [x-force y-force]
% FORCE(3) = struct('shape',forceshape,'values',values, ...
%     'nodes',nodes_with_force,'coords',[]);

% define a structure, that defines the properties of the applied forces
% first set of tractions
FORCE(1) = struct('shape',[],'values',[],'nodes',[],'coords',[]);
FORCE(1).nodes = [4 10 9 3];
FORCE(1).values = [0 0; 0 -1];        % [x-force y-force]
FORCE(1).shape = 'linear';

% nodes_with_force = [2 3 5 6 7];
% values = [10 -10];        % [x-force y-force]
% forceshape = 'constant';
% 
% FORCE(1) = struct('shape',forceshape,'values',values,'NBCnodes',nodes_with_force);

% % second set of tractions
% nodes_with_force = [4 14 15 16 3];
% values = [0 10];        % [x-force y-force]
% forceshape = 'constant';
% 
% FORCE(2) = struct('shape',forceshape,'values',values,'NBCnodes',nodes_with_force);
% 
% % third set of tractions
% nodes_with_force = [1 10 9 8 2];
% values = [0 -10];        % [x-force y-force]
% forceshape = 'constant';
% 
% FORCE(3) = struct('shape',forceshape,'values',values,'NBCnodes',nodes_with_force);
