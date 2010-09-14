% multi1_DBC.m
%
% Dirichlet boundary conditions for a patch test problem to check the multigrain
% example with 23 grains.

% Left hand side zero constraints
dispbc(1,1) = 1.0;
dispbc(1,2) = 1.0;
dispbc(1,3) = 1.0;
dispbc(1,4) = 1.0;
dispbc(1,5) = 1.0;
dispbc(1,6) = 1.0;
dispbc(1,7) = 1.0;
dispbc(1,8) = 1.0;
dispbc(1,9) = 1.0;
dispbc(1,10) = 1.0;
dispbc(1,11) = 1.0;

dispbc(2,1) = 1.0;
dispbc(2,2) = 1.0;
dispbc(2,3) = 1.0;
dispbc(2,4) = 1.0;
dispbc(2,5) = 1.0;
dispbc(2,6) = 1.0;
dispbc(2,7) = 1.0;
dispbc(2,8) = 1.0;
dispbc(2,9) = 1.0;
dispbc(2,10) = 1.0;
dispbc(2,11) = 1.0;

for i=441:451
  dispbc(1,i) = 1;
  ubar(1,i) = 1.;
end;


% % define enriched surfaces for bcs 
% bc_enr = 1;     % switch to show, that enriched DOFs are under BCs
% 
% 
% num_enr_surf = 6;
% 
% enr_surfs(1).nodes  = [2 3];
% enr_surfs(1).xsi    = [-1.0 1.0];   % This is in 1D xsi coords 
%                                     % related to the segment nodes
%                                     
% enr_surfs(1).coords = [0.0 8/5;     % There are the real coordinates of
%                        0.0 6/5];    % the subsegment
%                    
% enr_surfs(1).grain  = 1;
% enr_surfs(1).xy = [1 1];
% 
% enr_surfs(2).nodes  = [3 4];
% enr_surfs(2).xsi    = [-1.0 -0.5];
% enr_surfs(2).coords = [0.0 6/5;
%                        0.0 1.1];                  
% enr_surfs(2).grain  = 1;
% enr_surfs(2).xy = [1 1];
% 
% enr_surfs(3).nodes  = [3 4];
% enr_surfs(3).xsi    = [-0.5 1.0];
% enr_surfs(3).coords = [0.0 1.1;
%                        0.0 4/5];
% enr_surfs(3).grain  = 2;
% enr_surfs(3).xy = [1 1];
% 
% enr_surfs(4).nodes  = [7 8];
% enr_surfs(4).xsi    = [-1.0 1.0];
% enr_surfs(4).coords = [0.0 -2/5;
%                        0.0 -4/5];
% enr_surfs(4).grain  = 2;
% enr_surfs(4).xy = [1 1];
% 
% enr_surfs(5).nodes  = [8 9]
% enr_surfs(5).xsi    = [-1.0 -0.5];
% enr_surfs(5).coords = [0.0 -4/5;
%                        0.0 -0.9];
% enr_surfs(5).grain  = 2;
% enr_surfs(5).xy = [1 1];
% 
% enr_surfs(6).nodes  = [8 9]
% enr_surfs(6).xsi    = [-0.5 1.0];
% enr_surfs(6).coords = [0.0 -0.9;
%                        0.0 -6/5];
% enr_surfs(6).grain  = 3;
% enr_surfs(6).xy = [1 1];


