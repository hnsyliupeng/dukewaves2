% plotlagrangemultipliers.m
%
% This routine plots the values of lagrange multipliers over the interface.
% For each sliding case 'IFsliding_switch' there is an own routine, because
% the ordering and interpretation of the lagrange multipliers differs from
% case to case.
%

% Author: Matthias Mayr (04/2010)
    
% initialize some variables
xcoord = [];        % x-coords of all interface nodes
ycoord = [];        % y-coords of all interface nodes

% get index arrays of cut elements, sorted by number of cuts
cutelements1 = find(cutlist == 0);  % all uncut elements
cutelements2 = find(cutlist == 2);  % 1 interface in element
cutelements3 = find(cutlist == 3);  % two or three interfaces
cutelements = find(cutlist ~= 0);   % all cut elements        
        
% get max and min lagrange multipllier to scale the color range
maxlambda = max(lagmult);
minlambda = min(lagmult);

% scale the colormap
V = [minlambda maxlambda];
colmap = colormap(jet(100));
caxis(V);
colorbar;

% depending on the sliding case, different plots are made
switch IFsliding_switch  % different plot routines for sliding / no sliding
    case 0                      % no sliding at all (fully constrained)
        % create a new figure with three subplots
        figure(1);      
        subplot(311)    % absolute values of lagrange multipliers
        subplot(312);   % lagrange multipliers in normal direction
        subplot(313);   % lagrange multipliers in tangential direction
        hold on;
        
        % Plot interfaces without triple junctions first        
        for i=1:length(cutelements)
            for j=1:numgrain                % j = grain ID
                if seg_cut_info(j,i).nb_int == 2    % only elemetens with 
                                                    % one interface
                    % get current element ID 'eleID'
                    eleID = seg_cut_info(j,i).elemno;
                    
                    % get x- and y-coords of intersection points
                    xcoord = seg_cut_info(j,i).xint(:,1);
                    ycoord = seg_cut_info(j,i).xint(:,2);
                    
                    % get row in 'lagmult', in which lagrange multipliers for this
                    % interface start
                    lagmult_row = ...
                        (lag_surf(find(lag_surf(:,2)==eleID),1) * 2) - 1;
                    
                    % get x- and y-values of lagrange multiplier in element 'i'
                    lag_x = lagmult(lagmult_row);
                    lag_y = lagmult(lagmult_row + 1);

                    % compute absolute value
                    lag_val = sqrt(lag_x^2 + lag_y^2);
                    
                    % get tangential vector
                    tangvec = [xcoord(2)-xcoord(1);ycoord(2)-ycoord(1)];
                    tangvec = tangvec/norm(tangvec);
                    
                    % get normal vector
                    normvec = [-tangvec(2);tangvec(1)];
                    
                    % Now, there are three plots
                    % (1) Absolute values of lagrange multipliers
                    subplot(311);
                    % define line color (scaled to lagrange multiplier
                    colindex = ...
                        ceil((lag_val-minlambda)/(maxlambda-minlambda)*100);
                    % limit 'colindex' to values 1 <= colindex <= 100
                    if colindex < 1         
                        colindex = 1; 
                    elseif colindex > 100
                        colindex = 100;
                    end;
                    linecolor = colmap(colindex,:);

                    % plot interface
                    line(xcoord,ycoord,'Color',linecolor,'LineWidth',3);
                    
                    % (2) lagrange multipliers in normal direction
                    subplot(312);
                    % compute lagrange multiplier in normal direction
                    lag_normal = [lag_x lag_y] * normvec;
                    
                    % define line color (scaled to lagrange multiplier
                    colindex = ...
                        ceil((lag_normal-minlambda)/(maxlambda-minlambda)*100);
                    % limit 'colindex' to values 1 <= colindex <= 100
                    if colindex < 1         
                        colindex = 1; 
                    elseif colindex > 100
                        colindex = 100;
                    end;
                    linecolor = colmap(colindex,:);

                    % plot interface
                    line(xcoord,ycoord,'Color',linecolor,'LineWidth',3);
                    
                    % (3) lagrange multipliers in tangential direction
                    subplot(313);
                    % compute lagrange multiplier in tangential direction
                    lag_tang = [lag_x lag_y] * tangvec;
                    
                    % define line color (scaled to lagrange multiplier
                    colindex = ...
                        ceil((lag_tang-minlambda)/(maxlambda-minlambda)*100);
                    % limit 'colindex' to values 1 <= colindex <= 100
                    if colindex < 1         
                        colindex = 1; 
                    elseif colindex > 100
                        colindex = 100;
                    end;
                    linecolor = colmap(colindex,:);

                    % plot interface
                    line(xcoord,ycoord,'Color',linecolor,'LineWidth',3);
                elseif seg_cut_info(j,i).nb_int > 2
                    error('Can´t plot lagrange multipliers for triple junctions, yet.');
                end;
            end;
        end;
                
    case 1                      % frictionless sliding
        % For frictionless sliding, only one equation per constraint will
        % be added to 'bigk'. So, some indices change. 
        % Due to the frictionless sliding, there is no tangential lagrange
        % multiplier. So, only one plot will be generated.
        
        % create a new figure (no subplots due to frictionless sliding)
        figure(1);      
        hold on;
        
        % Plot interfaces without triple junctions first        
        for i=1:length(cutelements)
            for j=1:numgrain                % j = grain ID
                if seg_cut_info(j,i).nb_int == 2    % only elemetens with 
                                                    % one interface
                    % get current element ID 'eleID'
                    eleID = seg_cut_info(j,i).elemno;
                    
                    % get x- and y-coords of intersection points
                    xcoord = seg_cut_info(j,i).xint(:,1);
                    ycoord = seg_cut_info(j,i).xint(:,2);
                    
                    % get row in 'lagmult', in which lagrange multipliers for this
                    % interface start
                    lagmult_row = lag_surf(find(lag_surf(:,2)==eleID),1);
                    
                    % get x- and y-values of lagrange multiplier in element 'i'
                    lag_normal = lagmult(lagmult_row);
                    
                     % define line color (scaled to lagrange multiplier
                    colindex = ...
                        ceil((lag_normal-minlambda)/(maxlambda-minlambda)*100);
                    % limit 'colindex' to values 1 <= colindex <= 100
                    if colindex < 1         
                        colindex = 1; 
                    elseif colindex > 100
                        colindex = 100;
                    end;
                    linecolor = colmap(colindex,:);

                    % plot interface
                    line(xcoord,ycoord,'Color',linecolor,'LineWidth',3);
                elseif seg_cut_info(j,i).nb_int > 2
                    error('Can´t plot lagrange multipliers for triple junctions, yet.');
                end;
            end;
        end;            
    otherwise
        error('MATLAB:XFEM:UnvalidID',...
            'Unvalid sliding ID. Choose a valid ID or add an addition case to switch-case-structure.');
end;

% print max and min lagrange multiplier in console
text = ['Max. lambda:   ' num2str(maxlambda) '  Min. lambda:    '...
    num2str(minlambda) '  Range:  ' num2str(maxlambda-minlambda)];
disp(text);