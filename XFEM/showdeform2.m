% showdeform2.m
%
% plot original and deformed mesh
%

disp('showdeform2 ...');

x_def = [];
y_def = [];

% create a new figure
figure(10);
hold on;
set(10,'Name','showdeform2');

% plot initial mesh
for e=1:numele
   plot(x(node(1:3,e)),y(node(1:3,e)),'Color',[0.8 0.8 0.8],...
       'LineWidth',0.1)
   plot([x(node(3,e)) x(node(1,e))],[y(node(3,e)) y(node(1,e))],...
       'Color',[0.8 0.8 0.8],'LineWidth',1.0)
end;


% get coordinates of deformed mesh
for i = 1:numnod
    x_def(i) = x(i) + dis(2*i-1);
    y_def(i) = y(i) + dis(2*i); 
end

% plot deformed mesh
for e = 1:numele
   plot(x_def(node(1:3,e)),y_def(node(1:3,e)),'LineWidth',1);
   hold on
   plot([x_def(node(3,e)) x_def(node(1,e))], ...
       [y_def(node(3,e)) y_def(node(1,e))],'LineWidth',1);
   hold on
end

%h = plot(x_def(1:numnod),y_def(1:numnod));
% set(h,'LineStyle','none');
% set(h,'Marker','o');
% set(h,'MarkerEdgeColor','Red');
% set(h,'MarkerFaceColor','Red');
% set(h,'MarkerSize',4)
% axis([-1 24 -7 7]);

% scale axes of plotted mesh
ax_x = (max(x) - min(x))/max(x);
ax_y = (max(y) - min(y))/max(y);
axis([min(x)-ax_x max(x)+ax_x min(y)-ax_y max(y)+ax_y]);
clear ax_x ax_y;


% hold on
% 
% g = plot(x_def(numnod_1+1:numnod),y_def(numnod_1+1:numnod));
% set(g,'LineStyle','none');
% set(g,'Marker','o');
% set(g,'MarkerEdgeColor','Green');
% set(g,'MarkerFaceColor','Green');
% set(g,'MarkerSize',3)
% axis([-10 60 -30 30]);
