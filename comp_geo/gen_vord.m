
function[vx, vy] = genvord(x)

%clear
%close all


[ vx, vy ] = voronoi_new ( x(:,1), x(:,2) );
%scatter ( x(:,1), x(:,2), 'b', 'filled' );




%k = dsearch ( p(:,1), p(:,2), t, xs, ys );

  %axis ( [ -1, 2, -1, 2 ] )
  %title ( 'Diagram using "Voronoi\_New"' );
  %hold off

  
%[v,c]=voronoin(x); 
%for i = 1:length(c) 
%if all(c{i}~=1)   % If at least one of the indices is 1, 
                  % then it is an open region and we can't 
                  % patch that.
%                  i
%   patch(v(c{i},1),v(c{i},2),i); % use color i.
%end
%end
%axis equal
%figure(1)
%hold on
%for i=1:length(c)
%    plot(x(i,1),x(i,2),'o')
%end
%figure(2)
%voronoi(x(:,1),x(:,2) )
