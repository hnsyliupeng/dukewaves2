clear
close all

p =[      0    2.0000;
          0    0.0000;
          0   -2.0000;
     2.0000    1.1000;
     2.0000   -1.0000;
     4.0000    2.0000;
     4.0000   -0.1000;
     4.0000   -2.0000;
     6.0000    1.0000;
     6.0000   -1.1000;
     8.0000    2.0000;
     8.0000    0.2000;
     8.0000   -2.0000;
    10.0000    1.0000;
    10.0000   -1.1000;
    12.0000    2.0000;
    12.0000   -0.1000;
    12.0000   -2.0000;
    14.0000    1.0000;
    14.0000   -1.1000;
    16.0000    2.0000;
    16.0000    0.0000;
    16.0000   -2.0000];

%load multi1.mat vx vy       % load polycyrstalline grain mesh
[vx, vy] = gen_vord(p);

% Run through points in vx and vy, which are end points of line segments,
% and clip them to the boundaries of the domain

len = size(vx,2);

for i = 1:2
    for j = 1:len
        if vx(i,j) < 0.0
            
            if i == 1
                a = 2;
            else
                a = 1;
            end
            
            x = 0.0;
            x1 = vx(i,j);
            x2 = vx(a,j);
            y1 = vy(i,j);
            y2 = vy(a,j);
            
            % interpolate the proper y value
            y = (x - x1)*(y2 - y1)/(x2 - x1) + y1;
            
            vx(i,j) = x;
            vy(i,j) = y;
                       
        elseif vx(i,j) > 16.0
                                   
            if i == 1
                a = 2;
            else
                a = 1;
            end
            
            x = 16.0;
            x1 = vx(i,j);
            x2 = vx(a,j);
            y1 = vy(i,j);
            y2 = vy(a,j);
            
            % interpolate the proper y value
            y = (x - x1)*(y2 - y1)/(x2 - x1) + y1;
            
            vx(i,j) = x;
            vy(i,j) = y;
                       
        end
        
        
        if vy(i,j) < -2.0
                        
            if i == 1
                a = 2;
            else
                a = 1;
            end
            
            y = -2.0;
            x1 = vx(i,j);
            x2 = vx(a,j);
            y1 = vy(i,j);
            y2 = vy(a,j);
            
            x = (y - y1)*(x2 - x1)/(y2 - y1) + x1;
     
            vy(i,j) = y;
            vx(i,j) = x;
            
        elseif vy(i,j) > 2.0
                                    
            if i == 1
                a = 2;
            else
                a = 1;
            end
            
            y = 2.0;
            x1 = vx(i,j);
            x2 = vx(a,j);
            y1 = vy(i,j);
            y2 = vy(a,j);
            
            x = (y - y1)*(x2 - x1)/(y2 - y1) + x1;
     
            vy(i,j) = y;
            vx(i,j) = x;

        end
    end
end

for i = len:-1:1        %Get rid of some points
    if abs(vx(1,i) - 2.0) < 0.21
        vx(:,i) = [];
        vy(:,i) = [];
    elseif abs(vx(1,i) - 4.0) < 0.21
        vx(:,i) = [];
        vy(:,i) = [];
    elseif abs(vx(1,i) - 6.0) < 0.21
        vx(:,i) = [];
        vy(:,i) = [];
    elseif abs(vx(1,i) - 8.0) < 0.21
        vx(:,i) = [];
        vy(:,i) = [];
    elseif abs(vx(1,i) - 10.0) < 0.21
        vx(:,i) = [];
        vy(:,i) = [];
    elseif abs(vx(1,i) - 12.0) < 0.21
        vx(:,i) = [];
        vy(:,i) = [];
    elseif abs(vx(1,i) - 14.0) < 0.21
        vx(:,i) = [];
        vy(:,i) = [];
    elseif (abs(vx(1,i) - 16.0) < 0.1) && (abs(vy(1,i) - 0.0) < 0.21)
        vx(:,i) = [];
        vy(:,i) = [];
    elseif (abs(vx(1,i) - 0.0) < 0.1) && (abs(vy(1,i) - 0.0) < 0.21)
        vx(:,i) = [];
        vy(:,i) = [];
    end
end

        
x = [];
y = [];
x = [p(:,1)' vx(1,:) vx(2,:)];
y = [p(:,2)' vy(1,:) vy(2,:)];

% Seek and destroy duplicate points
% This is really a really inefficient routine, but... whatever

len = size(x,2);

for i = len:-1:1
    len2 = size(x,2);
    for j = len2:-1:1
        if (i ~= j) &&...
           ((abs(x(i) - x(j)) < 0.000001) &&...
            (abs(y(i) - y(j)) < 0.000001))
           
                x(i) = [];
                y(i) = [];
                break
        end
    end
end               

node = delaunay(x,y)'

numele = size(node,2);

% figure(1)
% hold on
% for e=1:numele
%    plot(x(node(1:3,e)),y(node(1:3,e)))
%    plot([x(node(3,e)) x(node(1,e))],[y(node(3,e)) y(node(1,e))])
% end
% 
% plot(vx,vy,'r','lineWidth',3)

% Associate an element with a grain and a grain-key-node

for i = 1:numele
    tag_element(i) = min(node(:,i));
end

% Associate all nodes with a grain or with grain 0 if they are on the edge.

tag_node = zeros(size(x,2),1);
tag_node(1:23) = [1:23];

% Loop through the elements, and create new nodes at all segment midpoints.
%  If the segment is associated with a grain-key-node, associate the new
%  node with that grain as well.  If not, tag it with grain zero.

cnt = size(x,2);                % New node number

for i = 1:numele                % Loop through elements
    
    key = tag_element(i);
    
    for j = 1:3                 % Loop through edges
        
        cnt = cnt + 1;
        
        if j == 1                % First edge
        
            n1 = node(1,i);
            n2 = node(2,i);
            x1 = x(n1);
            x2 = x(n2);
            y1 = y(n1);
            y2 = y(n2);
            x(cnt) = (x1 + x2)/2;
            y(cnt) = (y1 + y2)/2;
            if (n1 == key) || (n2 == key)
                tag_node(cnt) = key;
            else
                tag_node(cnt) = 0;
            end
            
        elseif j == 2                % Second edge
        
            n1 = node(2,i);
            n2 = node(3,i);
            x1 = x(n1);
            x2 = x(n2);
            y1 = y(n1);
            y2 = y(n2);
            x(cnt) = (x1 + x2)/2;
            y(cnt) = (y1 + y2)/2;
            if (n1 == key) || (n2 == key)
                tag_node(cnt) = key;
            else
                tag_node(cnt) = 0;
            end

        elseif j == 3                % Second edge
        
            n1 = node(3,i);
            n2 = node(1,i);
            x1 = x(n1);
            x2 = x(n2);
            y1 = y(n1);
            y2 = y(n2);
            x(cnt) = (x1 + x2)/2;
            y(cnt) = (y1 + y2)/2;
            if (n1 == key) || (n2 == key)
                tag_node(cnt) = key;
            else
                tag_node(cnt) = 0;
            end            
        end
    end
end
            
% Seek and destroy duplicate points
% This is really a really inefficient routine, but... whatever

len = size(x,2);

for i = len:-1:1
    len2 = size(x,2);
    for j = len2:-1:1
        if (i ~= j) &&...
           ((abs(x(i) - x(j)) < 0.00001) &&...
            (abs(y(i) - y(j)) < 0.00001))
           
                x(i) = [];
                y(i) = [];
                tag_node(i) = [];
                break
        end
    end
end               
    
figure(1)
hold on
for e=1:numele
   plot(x(node(1:3,e)),y(node(1:3,e)))
   plot([x(node(3,e)) x(node(1,e))],[y(node(3,e)) y(node(1,e))])
end

plot(vx,vy,'r','lineWidth',3)
       
node2 = delaunay(x,y)';

numele2 = size(node2,2);
        
figure(2)
hold on
for e=1:numele2
   plot(x(node2(1:3,e)),y(node2(1:3,e)))
   plot([x(node2(3,e)) x(node2(1,e))],[y(node2(3,e)) y(node2(1,e))])
end

plot(vx,vy,'r','lineWidth',3)

% Assign grains to elements

tag_element = zeros(numele2,1);

for i = 1:numele2
    for j = 1:3
        if tag_node(node2(j,i)) ~= 0
            i;
            tag_element(i) = tag_node(node2(j,i));
            break
        end
    end
end

figure(3)
hold on
for e=1:numele2
   patch(x(node2(1:3,e)),y(node2(1:3,e)),tag_element(e))
%   plot([x(node2(3,e)) x(node2(1,e))],[y(node2(3,e)) y(node2(1,e))])
end
        
plot(vx,vy,'r','lineWidth',3)     

numele = numele2;
node = node2;
numnod = size(x,2);

%save cmulti1.mat x y node numele numnod tag_element
        
        
        
        
        
        
        
        
        
        
        
        