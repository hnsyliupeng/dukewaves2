% plotslidestateevolutionNewton.m
%
% Plots the evolution of the stick-slip-area during a Newton iteration.
%

% Author: Matthias Mayr (07/2010)

% create a new figure
figure(40);
hold on;

if iter == 1
  clf;  % clear current figure
%   close 40
end;

% plot interfaces
% loop over all interfaces
for i=1:size(seg_cut_info,1)
  % loop over all elements
  for e=1:size(seg_cut_info,2)
    % only cut elements
    if seg_cut_info(i,e).elemno ~= -1
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

      % set color depending on current slidestate
      stylecell = {'b','r'};

      if any(seg_cut_info(i,e).f_trial > 0)
        index = 2;
      else
        index = 1;
      end;
      
%       % get 'index' for stylecell, depending on the 'slidestate'
%       index = seg_cut_info(i,e).slidestate + 1;

      % check, if it is a horizontal or a vertical interface
      if abs(xcoord(1) - xcoord(2)) < abs(ycoord(1)-ycoord(2))
        % vertical interface
        plot([iter iter],ycoord,stylecell{index},'LineWidth',3.0);   
        xlabel('normalized pseudo-time');
        ylabel('iteration step');
      else
        % horizontal interface
        plot(xcoord,[iter iter],stylecell{index},'LineWidth',3.0);       
        xlabel('x-coordinate');
        ylabel('iteration step');
      end;
    end;
  end;
end;

% edit figure
% legend('stick','slip');
title('blue = stick, red = slip');
hold off;


