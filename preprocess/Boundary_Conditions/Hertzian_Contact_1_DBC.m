% DBCs

% fixed nodes
% fixed in x-direction
for i=32:66
    dispbc(1,i) = 1;
    ubar(1,i) = 0;
end;

dispbc(1,1) = 1;
ubar(1,1) = 0;

dispbc(1,5) = 1;
ubar(1,5) = 0;

% fixed in y-direction
dispbc(2,1) = 1;
ubar(2,1) = 0;

for i=67:81
    dispbc(2,i) = 1;
    ubar(2,i) = 0;
end;

dispbc(2,2) = 1;
ubar(2,2) = 0;

% fix enriched DOFs
bc_enr = 1;     % switch to show, that enriched DOFs are under BCs

num_enr_surf = 2;

enr_surfs(1).nodes  = [57 58];
enr_surfs(1).xsi    = [-1.0 -0.582956228421145];   % This is in 1D xsi coords 
                                    % related to the segment nodes
enr_surfs(1).coords = [0.0 3.063888888896114;     % There are the real coordinates of
                       0.0 3.0000001];    % the subsegment
enr_surfs(1).grain  = 1;
enr_surfs(1).xy = [1 0];

enr_surfs(2).nodes  = [57 58];
enr_surfs(2).xsi    = [-0.582956228421145 1.0];   % This is in 1D xsi coords 
                                    % related to the segment nodes
enr_surfs(2).coords = [0.0 3.0000001;     % There are the real coordinates of
                       0.0 2.757500000006862];    % the subsegment
enr_surfs(2).grain  = 1;
enr_surfs(2).xy = [1 0];

