% minstabiparameter.m
%
% CALL: minstabiparameter(xcoords,ycoords,seg_cut_info,youngs, ...
%   endpoints,nodegrainmap,IFintegral)
%
% This routine is needed for Nitsche's method only. It estimates a minimum
% values for the stabilization parameter "alpha = penalty" in order to
% guarantee the positive definiteness of the system. For linear triangular
% elements, an estimate is given in "Dolbow, John and Harari, Isaac: An 
% efficient finite element method for embedded interface problems. Int. J.
% Numer. Meth. Engng. 2009; 78(2): 229-252"
%
% Input parameters:
%   xcoords         x-coordinates of the element's nodes
%   ycoords         y-coordinates of the element's nodes
%   seg_cut_info    some data about the interface in that element
%   youngs          Young's moduli of the grains in that element
%   endpoints       endpoints, defining the interfaces cutting this element
%   nodegrainmap    stores, in which grain each node resides
%   IFintegral      variant of stabilization (influences the scaling of
%                   'h')
%
% Returned variables:
%   alpha_min       minimal value for stabilization parameter to guarantee
%                   positive definiteness

% Author: Matthias Mayr (07/2010)

function [alpha_min] = minstabiparameter(xcoords,ycoords,seg_cut_info, ...
  youngs,endpoints,nodegrainmap,IFintegral)

% initialize
area_1 = 0;
area_2 = 0;

% get intersection points of interface with element edges
intersection = seg_cut_info.xint;

% end points of intersection - direction doesn't matter - this is for the
% segment jacobian calculation
if all(size(intersection) == [2 2])     % no triple junction in the element
  % get the two intersection points
  p1 = intersection(1,:);
  p2 = intersection(2,:);
elseif all(size(intersection) == [1 2]) % triple junction in the element
  % get the two points, that define the subsegment
  p1 = intersection(1,:);

  % Second endpoint of segment is also end point of interface
  endpoint = endpoints(1,:);

  inside = polygon_contains_point_2d ( 3, [xcoords;ycoords], endpoint );

  if inside
    p2 = endpoint;
  else
    p2 = endpoints(2,:);
  end
end

% for i=1:3
%   v = [p1' p2' [xcoords(i);ycoords(i)]];
%   currentarea = polygon_area_2d_2(3,v);
%   
%   if nodegrainmap(i) == seg_cut_info.grains(1)
%     area_1 = area_1 + abs(currentarea);
%   else
%     area_2 = area_2 + abs(currentarea);
%   end;
% end;

% area of the entire element
elementarea = abs(det([[1 1 1]' xcoords' ycoords']) / 2);

% indices of nodes in each grains
nodesgrain1 = find(nodegrainmap == seg_cut_info.grains(1));
nodesgrain2 = find(nodegrainmap == seg_cut_info.grains(2));

% choose the grain with only one node in it and get the nodes coordinates
if nodesgrain1 == 1
  xnode = xcoords(nodesgrain1(1));
  ynode = ycoords(nodesgrain1(1));
else
  xnode = xcoords(nodesgrain2(1));
  ynode = ycoords(nodesgrain2(1));
end;

% add the intersection points to the coordinate arrays
xnode = [xnode p1(1) p2(1)];
ynode = [ynode p1(2) p2(2)];

% compute the area of the part of the cut element which is a triangle
littlearea = abs(triangle_area_2d([xnode; ynode]));

% choose the grain with only one node in it and compute the areas of the
% two parts of the element covering the two grains
if nodesgrain1 == 1
  area_1 = littlearea;
  area_2 = elementarea - littlearea;
else
  area_1 = elementarea - littlearea;
  area_2 = littlearea;
end;

% length of subsegment of interface
Lsubsegment = norm(p1 - p2);

% compute alpha_min similar to eq. (42) and (53) in the paper mentioned
% above, but depending on the scaling of alpha by 'h'
switch IFintegral 
  case 1
    % compute alpha_min using eq. (42) and (53) in the paper mentioned above
    alpha_min = Lsubsegment / 2 * (youngs(1) / area_1 + youngs(2) / area_2);
  case 2
    alpha_min = Lsubsegment^2 / 2 * (youngs(1) / area_1^2 + youngs(2) / area_2^2);
  otherwise
    error('MATLAB:XFEM:UnvalidID', ...
      'Unvalid ID for number of integrals in penalty term');
end;

end

