% plottraction_frictionless_sliding1_plasticity.m
%
% Plots the tangential traction for the example
% 'inp_frictionless_sliding1_72_18.m' when perfect plasticity is introduced
% at the interface.
%

% Author: Matthias Mayr (08/2010)

disp('plottraction ...');

%% INITIALIZE
maxnormal = -inf;
minnormal = +inf;

maxtangent = -inf;
mintangent = +inf;

plotdata = [];
% ----------------------------------------------------------------------- %
%% FIND MAXIMA AND MINIMA
% loop over all subsegments and find the maximal and minimal traction value
% in normal and tangential direction.

% loop over all interfaces 'i'
for i=1:size(seg_cut_info,1)
  % loop over cut elements 'e'
  for e=1:size(seg_cut_info,2)
    % only, if 'e' is cut by 'i'
    if seg_cut_info(i,e).elemno ~= -1
      % get maximal and minimal values of current subsegment
      currentmaxnormal = max(seg_cut_info(i,e).ntracconv);  % MAX of normal
      currentminnormal = min(seg_cut_info(i,e).ntracconv);  % MIN of normal
      currentmaxtangent = max(seg_cut_info(i,e).ttracconv); % MAX of tangential
      currentmintangent = min(seg_cut_info(i,e).ttracconv); % MIN of tangential
    
      % compare 'currentmaxnormal' to 'maxnormal'
      if currentmaxnormal > maxnormal
        maxnormal = currentmaxnormal;
      end;
      
      % compare 'currentminnormal' to 'minnormal'
      if currentminnormal < minnormal
        minnormal = currentminnormal;
      end;
      
      % compare 'currentmaxtangent' to 'maxtangent'
      if currentmaxtangent > maxtangent
        maxtangent = currentmaxtangent;
      end;
      
      % compare 'currentmintangent' to 'mintangent'
      if currentmintangent < mintangent
        mintangent = currentmintangent;
      end;
    end;
  end;
end;
% Now, 'maxnormal', 'minnormal', 'maxtangent' and 'mintangent' contain the
% maximal and minimal values of normal and tangential tractions in the
% entire discretization.
% ----------------------------------------------------------------------- %
%% PRINT MAX AND MIN VALUES INTO CONSOLE
text = sprintf('Max_normal:\t\t%f\tMin_normal:\t\t%f\nMax_tangential:\t%f\tMin_tangential:\t%f', ...
  maxnormal,minnormal,maxtangent,mintangent); 
disp(text);
% ----------------------------------------------------------------------- %
%% PLOT NORMAL AND TANGENTIAL TRACTION AT GAUSS POINTS
% loop over all interfaces 'i'
for i=1:size(seg_cut_info,1)
  % loop over cut elements 'e'
  for e=1:size(seg_cut_info,2)
    % only, if 'e' is cut by 'i'
    if seg_cut_info(i,e).elemno ~= -1
      % PREPARE FINDING THE GAUSS POINTS
      % end points of intersection - direction doesn't matter - this is for the
      % segment jacobian calculation
      intersection = seg_cut_info(i,e).xint;
      if all(size(intersection) == [2 2])     % no triple junction in the element
        % get the two intersection points
        p1 = intersection(1,:);
        p2 = intersection(2,:);
      elseif all(size(intersection) == [1 2]) % triple junction in the element
        % get the two points, that define the subsegment
        p1 = intersection(1,:);

        % Second endpoint of segment is also end point of interface
        endpoint = INTERFACE_MAP(i).endpoints(1,:);
        
        % get some element data
        eleID = seg_cut_info(i,e).elemno;
        elenodes = node(:,eleID);
        xcoords = x(elenodes);
        ycoords = y(elenodes);
        
        % check, if 'endpoint' resides in element
        inside = polygon_contains_point_2d ( 3, [xcoords;ycoords], endpoint );

        if inside
          p2 = endpoint;
        else
          p2 = INTERFACE_MAP(i).endpoints(2,:);
        end
      end

      % Gauss points on segments
      gauss = [-sqrt(3)/3 sqrt(3)/3];

      % LOOP OVER GAUSS POINTS
      for g = 1:2
        % Get real coordinates of gauss point 'g'
        xn = 0.5*(1-gauss(g))*p1(1)+0.5*(1+gauss(g))*p2(1);
        yn = 0.5*(1-gauss(g))*p1(2)+0.5*(1+gauss(g))*p2(2);  
        
        % get values of normal traction at gauss point
        ntrac = seg_cut_info(i,e).ntracconv(g);
        ttrac = seg_cut_info(i,e).ttracconv(g);
        
%         figure(65)                      % normal traction
%         hold on;
%         set(65,'Name','normal traction');
%         subplot(211);
%         title('horizontal interface');
%         xlabel('x');
%         ylabel('traction');
%         hold on;
%         plot(xn,ntrac,'.');             % horizontal interface
%         hold off;
%         subplot(212);
%         title('vertical interface');
%         xlabel('traction');
%         ylabel('y');
%         hold on;
%         plot(ntrac,yn,'.');             % vertical interface
%         hold off;
        
%         figure(2)                      % tangential traction
%         hold on;
%         set(66,'Name','tangential traction');
%         subplot(211);
%         title('horizontal interface');
%         xlabel('x');
%         ylabel('traction');
%         hold on;
%         plot(xn,ttrac,'.');             % horizontal interface
%         hold off;
%         subplot(212);
%         title('vertical interface');
%         xlabel('traction');
%         ylabel('y-coordinate');
%         hold on;
        plotdata = [plotdata;
                    yn -ntrac];
%         plot(ttrac,yn,'.r');             % vertical interface
%         hold off;
      end;
    end;
  end;
end;

% sort 'plotdata' with respect to the y-coordinate of the Gauss points
plotdata_sort = sortrows(plotdata);

figure(46);
xlabel('normal traction','FontSize',36);
ylabel('y-coordinate','FontSize',36);
hold on;
plot(plotdata_sort(:,2),plotdata_sort(:,1),'b-','LineWidth',6);
axis([-6 6 -2 2])
% ----------------------------------------------------------------------- %
%% plot analytical solution
figure(46)
hold on;
plot([-1 1],[-2 2],'r-','LineWidth',3);
legend('simulation','analytical','Location','NorthWest')
axis([-1.5 1.5 -2 2]);
% ----------------------------------------------------------------------- %
%% clear some temporary variables
clear eleID elenodes xcoords ycoords i j endpoint p1 p2 xn yn ntrac ...
  ttrac text maxnormal minnormal maxtangent mintangent currentmaxnormal ...
  currentminnormal currentmaxtangent currentmintangent inside gauss g ...
  intersection plotdata plotdata_sort;
% ----------------------------------------------------------------------- %