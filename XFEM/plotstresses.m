% plotstresses.m
%
% Plot the xx-, yy- and xy-stresses. Subelements are taken in consideration
% for these plots.
%

% Author: Matthias Mayr (04/2010)

%% configure plot routine
% If you want plot the element edges, then set 'plotedges' to 1, else to 0
plotedges = 0;

%% get maximum dimensions of domain
xmax = max(X);
xmin = min(X);
ymax = max(Y);
ymin = min(Y);

%% plot xx-stresses
% open a new figure
figure(1);
subplot(311);
% axis equal;
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
colorbar;

% first, plot all uncut elements
% loop over all uncut elements
for e = 1:numele
    
    minix = x(node(:,e));
    miniy = y(node(:,e));
    
    if cutlist(e) == 0
        % get grain ID of element
        grain = SUBELEMENT_GRAIN_MAP(e);

        % plot colored element
        if plotedges == 1
            patch(minix,miniy,stress(e,4,grain));
        else
            patch(minix,miniy,stress(e,4,grain),'EdgeColor','none');
        end;
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
    if plotedges == 1
        patch(minix,miniy,stress(e,4,grain));
    else
        patch(minix,miniy,stress(e,4,grain),'EdgeColor','none');
    end;
end;

hold off;

%% plot yy-stresses
% open a new figure
% figure(2);
subplot(312);
% axis equal;
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
colorbar;

% first, plot all uncut elements
% loop over all uncut elements
for e = 1:numele
    
    minix = x(node(:,e));
    miniy = y(node(:,e));
    
    if cutlist(e) == 0
        % get grain ID of element
        grain = SUBELEMENT_GRAIN_MAP(e);

        % plot colored element
        if plotedges == 1
            patch(minix,miniy,stress(e,5,grain));
        else
            patch(minix,miniy,stress(e,5,grain),'EdgeColor','none');
        end;
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
    if plotedges == 1
        patch(minix,miniy,stress(e,5,grain));
    else
        patch(minix,miniy,stress(e,5,grain),'EdgeColor','none');
    end;
end;

hold off;

%% plot xy-stresses
% open a new figure
% figure(3);
subplot(313);
% axis equal;
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
colorbar;

% first, plot all uncut elements
% loop over all uncut elements
for e = 1:numele
    
    minix = x(node(:,e));
    miniy = y(node(:,e));
    
    if cutlist(e) == 0
        % get grain ID of element
        grain = SUBELEMENT_GRAIN_MAP(e);

        % plot colored element
        if plotedges == 1
            patch(minix,miniy,stress(e,6,grain));
        else
            patch(minix,miniy,stress(e,6,grain),'EdgeColor','none');
        end;
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
    if plotedges == 1
        patch(minix,miniy,stress(e,6,grain));
    else
        patch(minix,miniy,stress(e,6,grain),'EdgeColor','none');
    end;
end;

hold off;

%colorbar('location','west')

% print the maximum stresses into the console
text_xx = ['Max. xx-stress:  ' num2str(maxstress_xx) ...
    '   Min. xx-stress:    ' num2str(minstress_xx)];
text_yy = ['Max. yy-stress:  ' num2str(maxstress_yy) ...
    '   Min. yy-stress:    ' num2str(minstress_yy)];
text_xy = ['Max. xy-stress:  ' num2str(maxstress_xy) ...
    '   Min. xy-stress:    ' num2str(minstress_xy)];
disp(text_xx);
disp(text_yy);
disp(text_xy);

% clear some temporary variables
clear e grain minix miniy sub_nodes sube CONN_rows CONN_cols V plotedges;

