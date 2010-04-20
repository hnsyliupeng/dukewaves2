% plotstresses.m
%
% Plot the xx-, yy- and xy-stresses. Subelements are taken in consideration
% for these plots.
%

% Author: Matthias Mayr (04/2010)

%% get maximum dimensions of domain
xmax = max(X);
xmin = min(X);
ymax = max(Y);
ymin = min(Y);

%% plot xx-stresses
% open a new figure
figure(1);
subplot(311);
hold on;
title('xx-stress');
axis([xmin xmax ymin ymax]);

% find maximum and minimum stress
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

% scale the colormap
V = [minstress_xx maxstress_xx];
colormap(jet(100));
caxis(V);

% first, plot all uncut elements
% loop over all uncut elements
for e = 1:numele
    
    minix = x(node(:,e));
    miniy = y(node(:,e));
    
    if cutlist(e) == 0
        % get grain ID of element
        grain = SUBELEMENT_GRAIN_MAP(e);

        % plot colored element
        patch(minix,miniy,stress(e,4,grain));
    else
        
    end;
end;

% now, plot intersected/cut elements
% get dimensions of CONN
[CONN_rows CONN_cols] = size(CONN); % number of additional subelements = 
                                    % CONN_cols - numele
CONN_rows;
% loop over intersected/cut elements
for sube = numele+1:CONN_cols
    sub_nodes = CONN(:,sube);
    minix = X(sub_nodes);
	miniy = Y(sub_nodes);
    
    % get grain ID of element
    grain = SUBELEMENT_GRAIN_MAP(sube);
    
    % skip elements, which do not belong to a grain
    if grain == -1; continue; end;
    
    % get ID of parent element
    e=PARENTELEM_INFO(sube);
        
    % plot colored element
    patch(minix,miniy,stress(e,4,grain));
end;

hold off;

%% plot yy-stresses
% open a new figure
%figure(2);
subplot(312);
hold on;
title('yy-stress');
axis([xmin xmax ymin ymax]);

% find maximum and minimum stress
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

% scale the colormap
V = [minstress_yy maxstress_yy];
colormap(jet(100));
caxis(V);

% first, plot all uncut elements
% loop over all uncut elements
for e = 1:numele
    
    minix = x(node(:,e));
    miniy = y(node(:,e));
    
    if cutlist(e) == 0
        % get grain ID of element
        grain = SUBELEMENT_GRAIN_MAP(e);

        % plot colored element
        patch(minix,miniy,stress(e,5,grain));
    end;
end;

% now, plot intersected/cut elements
% get dimensions of CONN
[CONN_rows CONN_cols] = size(CONN); % number of additional subelements = 
                                    % CONN_cols - numele
CONN_rows;
for sube = numele+1:CONN_cols
    sub_nodes = CONN(:,sube);
    minix = X(sub_nodes);
	miniy = Y(sub_nodes);
    
    % get grain ID of element
    grain = SUBELEMENT_GRAIN_MAP(sube);
    
    % skip elements, which do not belong to a grain
    if grain == -1; continue; end;

    % get ID of parent element
    e=PARENTELEM_INFO(sube);
    
    % plot colored element
    patch(minix,miniy,stress(e,5,grain));
end;

hold off;

%% plot xy-stresses
% open a new figure
%figure(3);
subplot(313);
hold on;
title('xy-stress');
axis([xmin xmax ymin ymax]);

% find maximum and minimum stress
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

% scale the colormap
V = [minstress_xy maxstress_xy];
colormap(jet(100));
caxis(V);

% first, plot all uncut elements
% loop over all uncut elements
for e = 1:numele
    
    minix = x(node(:,e));
    miniy = y(node(:,e));
    
    if cutlist(e) == 0
        % get grain ID of element
        grain = SUBELEMENT_GRAIN_MAP(e);

        % plot colored element
        patch(minix,miniy,stress(e,6,grain));
    end;
end;

% now, plot intersected/cut elements
% get dimensions of CONN
[CONN_rows CONN_cols] = size(CONN); % number of additional subelements = 
                                    % CONN_cols - numele
for sube = numele+1:CONN_cols
    sub_nodes = CONN(:,sube);
    minix = X(sub_nodes);
	miniy = Y(sub_nodes);
    
    % get grain ID of element
    grain = SUBELEMENT_GRAIN_MAP(sube);

    % skip elements, which do not belong to a grain
    if grain == -1; continue; end;
    
    % get ID of parent element
    e=PARENTELEM_INFO(sube);
    
    % plot colored element
    patch(minix,miniy,stress(e,6,grain));
end;

hold off;

%colorbar('location','west')

