x_def = [];
y_def = [];

% % plot initial mesh
% for e=1:numele
%    plot(x(node(1:3,e)),y(node(1:3,e)),'Color',[0.5 0.5 0.5],...
%        'LineWidth',0.1)
%    plot([x(node(3,e)) x(node(1,e))],[y(node(3,e)) y(node(1,e))],...
%        'Color',[0.5 0.5 0.5],'LineWidth',1.0)
% end;


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
ax_x = (max(X) - min(X))/max(X);
ax_y = (max(Y) - min(Y))/max(Y);
axis([min(X)-ax_x max(X)+ax_x min(Y)-ax_y max(Y)+ax_y]);
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
