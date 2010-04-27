% plotlagrangemultipliers.m
%
% This routine plots the values of lagrange multipliers over the interface.
% For each sliding case 'IFsliding_switch' there is an own routine, because
% the ordering and interpretation of the lagrange multipliers differs from
% case to case.
%

% Author: Matthias Mayr (04/2010)
    
% get max and min lagrange multipllier to scale the color range
maxlambda = max(lagmult);
minlambda = min(lagmult);

if maxlambda == minlambda
    minlambda = 0.999999 * maxlambda;
end;

% scale the colormap
V = [minlambda maxlambda];
colmap = colormap(jet(100));
caxis(V);

% depending on the sliding case, different plots are made
switch IFsliding_switch  % different plot routines for sliding / no sliding
    case 0                      % no sliding at all (fully constrained)
        % create a new figure with three subplots
        figure(1);      
        subplot(311)    % absolute values of lagrange multipliers
        colorbar;
        axis equal;
        subplot(312);   % lagrange multipliers in normal direction
        colorbar;
        axis equal;
        subplot(313);   % lagrange multipliers in tangential direction
        colorbar;
        axis equal;
        hold on;
                
        % set figure window title
        set(1,'Name','Lagrange multipliers (no sliding)')
        
        %plot mesh into all subplots to have a better impression of 
        % the geometrie
        subplot(311)
        hold on
        title('absolute values');
        xlabel('x-coordinate');
        ylabel('y-coordinate');
        for e=1:numele
           plot(x(node(1:3,e)),y(node(1:3,e)),'Color',[0.9 0.9 0.9],...
               'LineWidth',0.1)
           plot([x(node(3,e)) x(node(1,e))],[y(node(3,e)) y(node(1,e))],...
               'Color',[0.9 0.9 0.9],'LineWidth',0.1)
        end;
        
        subplot(312)
        hold on
        title('normal direction');
        xlabel('x-coordinate');
        ylabel('y-coordinate');
        for e=1:numele
           plot(x(node(1:3,e)),y(node(1:3,e)),'Color',[0.9 0.9 0.9],...
               'LineWidth',0.1)
           plot([x(node(3,e)) x(node(1,e))],[y(node(3,e)) y(node(1,e))],...
               'Color',[0.9 0.9 0.9],'LineWidth',0.1)
        end;
        
        subplot(313)
        hold on
        title('tangential direction');
        xlabel('x-coordinate');
        ylabel('y-coordinate');
        for e=1:numele
           plot(x(node(1:3,e)),y(node(1:3,e)),'Color',[0.9 0.9 0.9],...
               'LineWidth',0.1)
           plot([x(node(3,e)) x(node(1,e))],[y(node(3,e)) y(node(1,e))],...
               'Color',[0.9 0.9 0.9],'LineWidth',0.1)
        end
        
        % loop over all interfaces
        for i = 1:size(seg_cut_info,1)      % every interface 'i'
            for e = 1:size(seg_cut_info,2)  % every cut element 'e' in 
                                            % interface 'i'
                if isempty(seg_cut_info(i,e).lagmult)==0    %only,if lagrange multiplier exists
                
                    % get 2 points, that determine the subsegment
                    if all(size(seg_cut_info(i,e).xint) == [2 2])       % no triple junction in 'e'
                        p1 = seg_cut_info(i,e).xint(1,:);
                        p2 = seg_cut_info(i,e).xint(2,:);
                    elseif all(size(seg_cut_info(i,e).xint) == [1 2])   % triple junction in 'e'
                        p1 = seg_cut_info(i,e).xint(1,:);

                        % Get coordinates of nodes of element
                        xep=zeros(1,3);
                        yep=zeros(1,3);
                        for m=1:3
                            jep = node(m,seg_cut_info(i,e).elemno); 
                            xep(m) = x(jep); 
                            yep(m) = y(jep);
                        end

                        % Second endpoint of segment is also end point of
                        % interface --> check, which one of the two endpoints
                        % of the interface lies in element 'e'

                        % get first endpoint
                        endpoint = INTERFACE_MAP(i).endpoints(1,:); 

                        % check, if it is in the element 'e'
                        inside = polygon_contains_point_2d ( 3, [xep;yep], endpoint );

                        if inside       % endpoint is in element
                            p2 = endpoint;
                        else            % endpoint is not in element
                            p2 = INTERFACE_MAP(i).endpoints(2,:);
                        end;
                    end;
                    % now, p1 and p2 are the two nodes, that determine the
                    % subsegment
                    xcoord = [p1(1) p2(1)];
                    ycoord = [p1(2) p2(2)];

                    % get lagrange multiplier for this subsegment (x and y)
                    lag = [seg_cut_info(i,e).lagmult(1) ...
                        seg_cut_info(i,e).lagmult(2)];

                    % compute absolute value, normal and tangential value
                    lag_val = norm(lag);
                    lag_normal = lag * seg_cut_info(i,e).normal;
                    lag_tang = lag * seg_cut_info(i,e).tangent;

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
        set(1,'Name','Lagrange multipliers (frictionless sliding');
        title('normal direction');
        xlabel('x-coordinate');
        ylabel('y-coordinate');
        colorbar;
        
        % plot mesh for a better impression of geometry
        hold on
        for e=1:numele
           plot(x(node(1:3,e)),y(node(1:3,e)),'Color',[0.9 0.9 0.9],...
               'LineWidth',0.1)
           plot([x(node(3,e)) x(node(1,e))],[y(node(3,e)) y(node(1,e))],...
               'Color',[0.9 0.9 0.9],'LineWidth',0.1)
        end;
            
        % loop over all interfaces
        for i = 1:size(seg_cut_info,1)      % every interface 'i'
            for e = 1:size(seg_cut_info,2)  % every cut element 'e' in 
                                            % interface 'i'
                if isempty(seg_cut_info(i,e).lagmult)==0    %only,if lagrange multiplier exists
                
                    % get 2 points, that determine the subsegment
                    if all(size(seg_cut_info(i,e).xint) == [2 2])       % no triple junction in 'e'
                        p1 = seg_cut_info(i,e).xint(1,:);
                        p2 = seg_cut_info(i,e).xint(2,:);
                    elseif all(size(seg_cut_info(i,e).xint) == [1 2])   % triple junction in 'e'
                        p1 = seg_cut_info(i,e).xint(1,:);

                        % Get coordinates of nodes of element
                        xep=zeros(1,3);
                        yep=zeros(1,3);
                        for m=1:3
                            jep = node(m,seg_cut_info(i,e).elemno); 
                            xep(m) = x(jep); 
                            yep(m) = y(jep);
                        end

                        % Second endpoint of segment is also end point of
                        % interface --> check, which one of the two endpoints
                        % of the interface lies in element 'e'

                        % get first endpoint
                        endpoint = INTERFACE_MAP(i).endpoints(1,:); 

                        % check, if it is in the element 'e'
                        inside = polygon_contains_point_2d ( 3, [xep;yep], endpoint );

                        if inside       % endpoint is in element
                            p2 = endpoint;
                        else            % endpoint is not in element
                            p2 = INTERFACE_MAP(i).endpoints(2,:);
                        end;
                    end;
                    % now, p1 and p2 are the two nodes, that determine the
                    % subsegment
                    xcoord = [p1(1) p2(1)];
                    ycoord = [p1(2) p2(2)];

                    % get lagrange multiplier for this subsegment (x and y)
                    lag = seg_cut_info(i,e).lagmult;

                    % compute absolute value, normal and tangential value
%                     lag_val = norm(lag);
                    lag_normal = lag * seg_cut_info(i,e).normal;
%                     lag_tang = lag * seg_cut_info(i,e).tangent;

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
                end;
            end;
        end;    
    otherwise
        error('MATLAB:XFEM:UnvalidID',...
            'Unvalid sliding ID. Choose a valid ID or add an addition case to switch-case-structure.');
end;

% figure(1);
% hold on;
% colorbar;
% hold off;

% print max and min lagrange multiplier in console
text = ['Max. lambda:   ' num2str(maxlambda) '  Min. lambda:    '...
    num2str(minlambda) '  Range:  ' num2str(maxlambda-minlambda)];
disp(text);

% clear some temporary variables
clear linecolor colmap xcoord ycoord lag_normal lag_tang lag_val lag_x ...
    lag_y text eleID lagmult_row colindex inside endpoint p1 p2 m jep ...
    xep yep i e lag;