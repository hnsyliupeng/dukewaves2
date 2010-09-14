% debugmaxstress.m
%

% Author: Matthias Mayr (09/2010)

%% find maximum and minimum stress
%% sigma_xx
maxstress_vec = max(stress(:,4,:));
minstress_vec = min(stress(:,4,:));

maxstress_xx = maxstress_vec(:,:,1);
minstress_xx = minstress_vec(:,:,1);
if maxngrains > 1
    for i=2:maxngrains
        if maxstress_vec(:,:,i) > maxstress_xx
            maxstress_xx = maxstress_vec(:,:,i);
        end;
        if minstress_vec(:,:,i) < minstress_xx
            minstress_xx = minstress_vec(:,:,i);
        end;
    end;
end;
clear maxstress_vec minstress_vec;
% ----------------------------------------------------------------------- %
%% sigma_yy
maxstress_vec = max(stress(:,5,:));
minstress_vec = min(stress(:,5,:));

maxstress_yy = maxstress_vec(:,:,1);
minstress_yy = minstress_vec(:,:,1);
if maxngrains > 1
    for i=2:maxngrains
        if maxstress_vec(:,:,i) > maxstress_yy
            maxstress_yy = maxstress_vec(:,:,i);
        end;
        if minstress_vec(:,:,i) < minstress_yy
            minstress_yy = minstress_vec(:,:,i);
        end;
    end;
end;
clear maxstress_vec minstress_vec;
% ----------------------------------------------------------------------- %
%% sigma_xy
maxstress_vec = max(stress(:,6,:));
minstress_vec = min(stress(:,6,:));

maxstress_xy = maxstress_vec(:,:,1);
minstress_xy = minstress_vec(:,:,1);
if maxngrains > 1
    for i=2:maxngrains
        if maxstress_vec(:,:,i) > maxstress_xy
            maxstress_xy = maxstress_vec(:,:,i);
        end;
        if minstress_vec(:,:,i) < minstress_xy
            minstress_xy = minstress_vec(:,:,i);
        end;
    end;
end;
clear maxstress_vec minstress_vec;
% ----------------------------------------------------------------------- %
%% von-Mises-stress
maxstress_vec = max(stressvonmises(:,4,:));
minstress_vec = min(stressvonmises(:,4,:));

maxstress_vonmises = maxstress_vec(:,:,1);
minstress_vonmises = minstress_vec(:,:,1);
if maxngrains > 1
    for i=2:maxngrains
        if maxstress_vec(:,:,i) > maxstress_vonmises
            maxstress_vonmises = maxstress_vec(:,:,i);
        end;
        if minstress_vec(:,:,i) < minstress_vonmises
            minstress_vonmises = minstress_vec(:,:,i);
        end;
    end;
end;
clear maxstress_vec minstress_vec;
% ----------------------------------------------------------------------- %
%% display
maxstress_xx
minstress_xx
maxstress_yy
minstress_yy
maxstress_xy
minstress_xy
maxstress_vonmises
minstress_vonmises