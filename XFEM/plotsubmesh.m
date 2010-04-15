close all
figure(1)
hold on
for e=1:num_sub_elems
   %plot(x(node(1:4,e))/4,y(node(1:4,e))/4+0.5)
   %plot([x(node(4,e))/4 x(node(1,e))/4],[y(node(4,e))/4+0.5 y(node(1,e))/4+0.5])
   plot(X(CONN(1:3,e)),Y(CONN(1:3,e)))
   plot([X(CONN(3,e)) X(CONN(1,e))],[Y(CONN(3,e)) Y(CONN(1,e))])
end