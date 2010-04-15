% boundary conditions for a patch test problem

ubar = zeros(2,numnod);
force = zeros(2,numnod);
dispbc = zeros(2,numnod);

dispbc(1,1) = 1.0;
dispbc(1,2) = 1.0;
dispbc(1,3) = 1.0;
%dispbc(1,4) = 1.0;
%dispbc(1,6) = 1.0;
dispbc(1,7) = 1.0;
dispbc(1,8) = 1.0;
dispbc(1,9) = 1.0;

dispbc(2,1) = 1.0;
dispbc(2,2) = 1.0;
dispbc(2,3) = 1.0;
%dispbc(2,4) = 1.0;
%dispbc(2,6) = 1.0;
dispbc(2,7) = 1.0;
dispbc(2,8) = 1.0;
dispbc(2,9) = 1.0;

ubar(1,1) = 0.0;
ubar(1,2) = 0.0;
ubar(1,3) = 0.0;
ubar(1,4) = 1.0;
ubar(1,6) = 1.0;
ubar(1,7) = 2.0;
ubar(1,8) = 2.0;
ubar(1,9) = 2.0;

ubar(2,1) = 0.0;
ubar(2,2) = 0.0;
ubar(2,3) = 0.0;
ubar(2,4) = 0.0;
ubar(2,6) = 0.0;
ubar(2,7) = 0.0;
ubar(2,8) = 0.0;
ubar(2,9) = 0.0;

% define enriched surfaces for bcs 

num_enr_surf = 6;

enr_surfs(1).nodes  = [6 9];
enr_surfs(1).xsi    = [-1.0 0.0];   % This is in 1D xsi coords 
                                    % related to the segment nodes
                                    
enr_surfs(1).coords = [1.0 0.0;     % There are the real coordinates of
                       1.5 0.0];    % the subsegment
                   
enr_surfs(1).grain  = 1;
enr_surfs(1).xy = [1 1];

enr_surfs(2).nodes  = [6 9];
enr_surfs(2).xsi    = [0.0 1.0];
enr_surfs(2).coords = [1.5 0.0;
                       2.0 0.0];                  
enr_surfs(2).grain  = 2;
enr_surfs(2).xy = [1 1];

enr_surfs(3).nodes  = [4 7];
enr_surfs(3).xsi    = [-1.0 0.0];
enr_surfs(3).coords = [1.0 2.0;
                       1.5 2.0];
enr_surfs(3).grain  = 1;
enr_surfs(3).xy = [1 1];

enr_surfs(4).nodes  = [4 7];
enr_surfs(4).xsi    = [0.0 1.0];
enr_surfs(4).coords = [1.5 2.0;
                       2.0 2.0];
enr_surfs(4).grain  = 2;
enr_surfs(4).xy = [1 1];

enr_surfs(5).nodes  = [1 4]
enr_surfs(5).xsi    = [-1 1];
enr_surfs(5).coords = [0.0 2.0;
                       1.0 2.0];
enr_surfs(5).grain  = 1;
enr_surfs(5).xy = [1 1];

enr_surfs(6).nodes  = [3 6]
enr_surfs(6).xsi    = [-1 1];
enr_surfs(6).coords = [0.0 0.0;
                       1.0 0.0];
enr_surfs(6).grain  = 1;
enr_surfs(6).xy = [1 1];
