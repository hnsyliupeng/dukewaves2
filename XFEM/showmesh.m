% showmesh.m
%
% plots the nodes.

disp('showmesh ...');

len = max(x);

h = plot(x(1:numnod),y(1:numnod));
set(h,'LineStyle','none');
set(h,'Marker','o');
set(h,'MarkerEdgeColor','Black');
set(h,'MarkerFaceColor','Black');
set(h,'MarkerSize',4)
axis([-1 len+2 -7 7]);
% 
% hold on
% 
% g = plot(x(numnod_1+1:numnod),y(numnod_1+1:numnod));
% set(g,'LineStyle','none');
% set(g,'Marker','o');
% set(g,'MarkerEdgeColor','Red');
% set(g,'MarkerFaceColor','Red');
% set(g,'MarkerSize',3)
% axis([-10 60 -30 30]);
