% DBCs

% fixed nodes
% fixed in x-direction
dispbc(1,1) = 1;
dispbc(1,2) = 1;

% fixed in y-direction
dispbc(2,1) = 1;
dispbc(2,2) = 1;




% fix enriched DOFs
bc_enr = 1;     % switch to show, that enriched DOFs are under BCs

num_enr_surf = 1;

enr_surfs(1).nodes  = [1 2];
enr_surfs(1).xsi    = [-1.0 1.0];   % This is in 1D xsi coords 
                                    % related to the segment nodes
enr_surfs(1).coords = [0.0 -2.0;     % There are the real coordinates of
                       0.0 2.0];    % the subsegment
enr_surfs(1).grain  = 1;
enr_surfs(1).xy = [1 1];

