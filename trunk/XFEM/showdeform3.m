% showdeform3.m
%
% Plot original mesh and deformed mesh in a way that the discontinuity can
% be seen [Routine adopted from 'showdeform2.m']
%

% Author: Matthias Mayr (08/2010)

disp('showdeform3 ...');

%% initialize
x_def = zeros(1,numnod);
y_def = zeros(1,numnod);

dispx = zeros(1,3);
dispy = zeros(1,3);

% create a new figure
figure(17);
hold on;
set(17,'Name','showdeform3');
clf;

% Each grain has its own color. This color is determined by the value of
% the grain number using a colormap. 
V = [1 maxngrains];
colormap(lines(maxngrains));
caxis(V);
% colorbar;
% ----------------------------------------------------------------------- %
%% plot initial mesh
% for e=1:numele
%    plot(x(node(1:3,e)),y(node(1:3,e)),'Color',[0.8 0.8 0.8],...
%        'LineWidth',0.1)
%    plot([x(node(3,e)) x(node(1,e))],[y(node(3,e)) y(node(1,e))],...
%        'Color',[0.8 0.8 0.8],'LineWidth',1.0)
% end;
% ----------------------------------------------------------------------- %
%% get coordinates of deformed mesh
for i = 1:numnod
    x_def(i) = x(i) + dis(2*i-1);
    y_def(i) = y(i) + dis(2*i); 
end
% ----------------------------------------------------------------------- %
%% plot deformed state
% plot only uncut elements
hold on;
for e = 1:numele
  if cutlist(e) == 0    % element is uncut
    % get grain ID of this subelement
    grain = find(elemgrainmap(e,:) == 1);
    
    % plot the element colored corresponding to the grain it belongs to
    patch(x_def(node(1:3,e)),y_def(node(1:3,e)),grain,'EdgeColor','none')
  else                  % element is cut
    % loop over all subelements of element 'e'
    for m = 1:SUBELEM_INFO(e).no_kids 
      % get global ID of subelement
      subele = SUBELEM_INFO(e).kids(m);

      % get grain ID of this subelement
      grain = SUBELEMENT_GRAIN_MAP(subele);

      % get nodes of subelement
      node_sub = CONN(:,subele);       

      % get coordinates of subelement nodes
      Xcoords = X(node_sub);
      Ycoords = Y(node_sub);

      % get data of parent element
      elenodes = node(:,e);
      xcoords = x(elenodes);
      ycoords = y(elenodes);

      % reset deformed coordinates of subelement nodes
      X_def = zeros(1,3);
      Y_def = zeros(1,3);

      if SUBELEMENT_GRAIN_MAP(subele) == 1
        % subelement is enriched

        % element displacement vector
        for b=1:3
            % base degrees of freedom
            b1 = id_eqns(elenodes(b),1);
            b2 = id_eqns(elenodes(b),2); 
            dispx(b) = totaldis(b1);
            dispy(b) = totaldis(b2);
        end

        % element displacement vector
        for b=1:3
            % extra degrees of freedom
            b1 = id_eqns(elenodes(b),3);
            b2 = id_eqns(elenodes(b),4); 
            dispx(b) = dispx(b) + totaldis(b1);
            dispy(b) = dispy(b) + totaldis(b2);
        end

        % get nodal displacements at parent nodes
        xd = dispx;
        yd = dispy;

        % compute area of parent element
        Area = det([[1 1 1]' xcoords' ycoords'])/2;

        % reset shape function matrix
        N = zeros(1,3);

        % interpolate the nodal displacements at the subelements nodes
        for n=1:3 % loop over the subelement's nodes
          xn = Xcoords(n);
          yn = Ycoords(n);
          
          N = zeros(1,3);
          
          % assemble 'N'
          for k=1:3
            % get area opposite of node of concern  
            xes = xcoords;
            yes = ycoords;
            xes(k) = xn;
            yes(k) = yn;
            Larea = det([[1 1 1]' xes' yes'])/2;
            
            N(k) = N(k) + Larea/Area;
          end;
          
          X_def(n) = Xcoords(n) + N * xd';
          Y_def(n) = Ycoords(n) + N * yd';
        end;

        patch(X_def,Y_def,grain,'EdgeColor','none')
      else
        % subelement is not enriched

        % get nodal displacements at parent nodes
        % element displacement vector
        for b=1:3
            % base degrees of freedom
            b1 = id_eqns(elenodes(b),1);
            b2 = id_eqns(elenodes(b),2); 
            dispx(b) = totaldis(b1);
            dispy(b) = totaldis(b2);
        end

        % get nodal displacements at parent nodes
        xd = dispx;
        yd = dispy;

        % compute area of parent element
        Area = det([[1 1 1]' xcoords' ycoords'])/2;

        % reset shape function matrix
        N = zeros(1,3);

        % interpolate the nodal displacements at the subelements nodes
        for n=1:3 % loop over the subelement's nodes
          xn = Xcoords(n);
          yn = Ycoords(n);
          
          N = zeros(1,3);
          
          % assemble 'N'
          for k=1:3
            % get area opposite of node of concern  
            xes = xcoords;
            yes = ycoords;
            xes(k) = xn;
            yes(k) = yn;
            Larea = det([[1 1 1]' xes' yes'])/2;
            
            N(k) = N(k) + Larea/Area;
          end;
          
          X_def(n) = Xcoords(n) + N * xd';
          Y_def(n) = Ycoords(n) + N * yd';
        end;
        
        patch(X_def,Y_def,grain,'EdgeColor','none')
      end;
    end;
  end;
end;


% clear some temporary variables
clear e i m subele node_sub dispx dispy b xd yd xe ye nodes grain eleID ...
  b1 b2 V n;
% ----------------------------------------------------------------------- %
%% scale axes of plotted mesh
ax_x = (max(x) - min(x))/max(x);
ax_y = (max(y) - min(y))/max(y);
axis([min(x)-ax_x max(x)+ax_x min(y)-ax_y max(y)+ax_y]);

% clear some temporary variables
clear ax_x ax_y;
% ----------------------------------------------------------------------- %
disp('finished');