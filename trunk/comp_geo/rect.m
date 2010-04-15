% Creates an unstructured tet mesh in a rectangular area
% Jessica Sanders
% 6/7/2006

beam_l = 20;
beam_h = 10;

ndiv_l = 6;
ndiv_h = 5;

% start wiith an evenly spaced mesh
for i = 1:(ndiv_l+1)
    for j = 1:(ndiv_h+1)
        x((ndiv_h+1)*(i-1)+j) = (beam_l/ndiv_l)*(i-1);
        y((ndiv_h+1)*(i-1)+j) = beam_h/2 - beam_h/ndiv_h*(j-1);
    end
end

variation_l = 0.2*beam_l/ndiv_l;
variation_h = 0.2*beam_h/ndiv_h;

for i = 2:(ndiv_l)
    for j = 2:(ndiv_h)
        x((ndiv_h+1)*(i-1)+j) = x((ndiv_h+1)*(i-1)+j) + randn*variation_l;
        y((ndiv_h+1)*(i-1)+j) = y((ndiv_h+1)*(i-1)+j) + randn*variation_h;
    end
end

node = (delaunay(x,y))';

numele = size(node,2);
numnod = size(x,2);
% 
% figure(1)
% hold on
% for e = 1:numele
%     plot(x(node(1:3,e)),y(node(1:3,e)));
%     plot([x(node(3,e)),x(node(1,e))],[y(node(3,e)),y(node(1,e))])
% end