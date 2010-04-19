x_def = [];
y_def = [];

for i = 1:numnod
    x_def(i) = x(i) + dis(2*i-1);
    y_def(i) = y(i) + dis(2*i); 
end

for e = 1:numele
   plot(x_def(node(1:3,e)),y_def(node(1:3,e)),'LineWidth',2);
   hold on
   plot([x_def(node(3,e)) x_def(node(1,e))],[y_def(node(3,e)) y_def(node(1,e))],'LineWidth',2);
   hold on
end

%h = plot(x_def(1:numnod),y_def(1:numnod));
% set(h,'LineStyle','none');
% set(h,'Marker','o');
% set(h,'MarkerEdgeColor','Red');
% set(h,'MarkerFaceColor','Red');
% set(h,'MarkerSize',4)
axis([-1 24 -7 7]);

% hold on
% 
% g = plot(x_def(numnod_1+1:numnod),y_def(numnod_1+1:numnod));
% set(g,'LineStyle','none');
% set(g,'Marker','o');
% set(g,'MarkerEdgeColor','Green');
% set(g,'MarkerFaceColor','Green');
% set(g,'MarkerSize',3)
% axis([-10 60 -30 30]);
