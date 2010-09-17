% moveintersectionpoints.m
% 
% CALL: moveintersectionpoints()
%
% Moves the intersection points of an interface with the element edges such
% that the smaller edge portion of a cut element has a length of at least 1
% percent of the edge length.
%
% Input parameters:
%   seg_cut_info      some data about interface subsegments
%   X                 x-coordinates of subelements nodes
%   Y                 y-coordinates of subelements nodes
%   x                 x-coordinates of elements nodes
%   y                 y-coordinates of elements nodes      
%   node              connectivity between nodes and elements
%   INTERFACE_MAP     infomration about interface geometry
%   tol               length tolerance to move intersection points
%
% Returned variables:
%   seg_cut_info      some data about interface subsegments (modified)
%   X                 x-coordinates of subelements nodes (modified)
%   Y                 y-coordinates of subelements nodes (modified)
%

% Author: Matthias Mayr (09/2010)

function [seg_cut_info X Y] = moveintersectionpoints(seg_cut_info, ...
  X,Y,x,y,node,INTERFACE_MAP,tol)

% loop over interfaces 'i'
for i=1:size(seg_cut_info,1)
  % loop over cut elements 'e'
  for e=1:size(seg_cut_info,2)
    % only if 'e' is cut by 'i'
    if seg_cut_info(i,e).elemno ~= -1
      % get some element data
      eleID = seg_cut_info(i,e).elemno;
      elenodes = node(:,eleID);
      xcoords = x(elenodes);
      ycoords = y(elenodes);
      
      endpoints = INTERFACE_MAP(i).endpoints;
      intersection = seg_cut_info(i,e).xint;
            
      % end points of intersection - direction doesn't matter - this is for the
      % segment jacobian calculation

      if all(size(intersection) == [2 2])
        p1 = intersection(1,:);
        p2 = intersection(2,:);
      elseif all(size(intersection) == [1 2])
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
      
      if seg_cut_info(i,e).nb_int == 1
        % This element contains a triple junction, so only 1 of the
        % intersection points lies on an element edge. The other one is
        % part of the triple junction and hence in the interior of the
        % element. Here, 'p1' is the intersection point. 'p2' does not have
        % to be considered here.
        
        % get edgeID of cut edge
        edgeID = seg_cut_info(i,e).edge_ids; % There should be only one
        
        switch edgeID
          case 1
            pe1 = [xcoords(1) ycoords(1)];
            pe2 = [xcoords(2) ycoords(2)];
          case 2
            pe1 = [xcoords(2) ycoords(2)];
            pe2 = [xcoords(3) ycoords(3)];
          case 3
            pe1 = [xcoords(3) ycoords(3)];
            pe2 = [xcoords(1) ycoords(1)];
          otherwise
            error('MATLAB:comp_geo','Unvalid edge ID');
        end;
        
        edge_length = norm(pe2 - pe1);
        
        if norm(p1 - pe1)/ edge_length < tol
          % 'p1' to close to 'pe1' --> move 'p1' towards 'pe2'
          p1 = (1 - tol) * pe1 + tol * pe2;          
          
        elseif norm(p1 - pe1)/ edge_length > (1 - tol)
          % 'p1' to close to 'pe2' --> move 'p1' towards 'pe1'
          p1 = tol * pe1 + (1 - tol) * pe2;
        else
          break
        end;
        
        % find subelement node, that should coincide with the intersection 
        % point
        diff_x = X - p1(1);
        diff_y = Y - p1(2);
        norm_vec = sqrt(diff_x.^2 + diff_y.^2);
        [val subnodeID] = min(norm_vec);
        
        if subnodeID < length(x)
            error('MATLAB:comp_geo','Wants to move a parent node');
        end;
        
        % update coordinates of subelement node, that should coincide with
        % the intersection point
        X(subnodeID) = p1(1);
        Y(subnodeID) = p1(2);
        
        % update intersection point in 'seg_cut_info'
        if inside
          seg_cut_info(i,e).xint(2,:) = p1;
        else
          seg_cut_info(i,e).xint(1,:) = p1;
        end;
      else
        % The element is cut only by one interface. Hence, there are two
        % intersection points which have to be checked. First, one has to
        % check, which intersection points belongs to which element edge.
        % Then, the same procedure as above can be applied.
        
        % assume, that 'p1' lies on the first cut element edge, check it
        % and correct, if necessary.
        p1_edge = seg_cut_info(i,e).edge_ids(1);
        switch p1_edge
          case 1
            pe1 = [xcoords(1) ycoords(1)];
            pe2 = [xcoords(2) ycoords(2)];
          case 2
            pe1 = [xcoords(2) ycoords(2)];
            pe2 = [xcoords(3) ycoords(3)];
          case 3
            pe1 = [xcoords(3) ycoords(3)];
            pe2 = [xcoords(1) ycoords(1)];
          otherwise
            error('MATLAB:comp_geo','Unvalid edge ID');
        end;
        dist_p1 = line_exp_point_dist_2d(pe1,pe2,p1);
        dist_p2 = line_exp_point_dist_2d(pe1,pe2,p2);
        
        % check assumption
        if dist_p1 < dist_p2
          % assumption right
          p2_edge = seg_cut_info(i,e).edge_ids(2);
        else
          % assumption wrong
          p1_edge = seg_cut_info(i,e).edge_ids(2);
          p2_edge = seg_cut_info(i,e).edge_ids(1);
        end;
        
        % build vector quantities to prepare a for-loop
        edges = [p1_edge p2_edge];
        int_points = [p1; p2];
        
        % loop over intersection points
        for g=1:2
          p1 = int_points(g,:);
          switch edges(g)
            case 1
              pe1 = [xcoords(1) ycoords(1)];
              pe2 = [xcoords(2) ycoords(2)];
            case 2
              pe1 = [xcoords(2) ycoords(2)];
              pe2 = [xcoords(3) ycoords(3)];
            case 3
              pe1 = [xcoords(3) ycoords(3)];
              pe2 = [xcoords(1) ycoords(1)];
            otherwise
              error('MATLAB:comp_geo','Unvalid edge ID');
          end;
          
          edge_length = norm(pe2 - pe1);
        
          if norm(p1 - pe1)/ edge_length < tol
            % 'p1' to close to 'pe1' --> move 'p1' towards 'pe2'
            p1 = (1 - tol) * pe1 + tol * pe2;
          elseif norm(p1 - pe1)/ edge_length > (1 - tol)
            % 'p1' to close to 'pe2' --> move 'p1' towards 'pe1'
            p1 = tol * pe1 + (1 - tol) * pe2;
          else
            break;
          end;

          % find subelement node, that should coincide with the intersection 
          % point
          diff_x = X - p1(1);
          diff_y = Y - p1(2);
          norm_vec = sqrt(diff_x.^2 + diff_y.^2);
          [val subnodeID] = min(norm_vec);

          if subnodeID < length(x)
            error('MATLAB:comp_geo','Wants to move a parent node');
          end;
          
          % update coordinates of subelement node, that should coincide with
          % the intersection point
          X(subnodeID) = p1(1);
          Y(subnodeID) = p1(2);   
          
          % update intersection points in 'seg_cut_info'
          seg_cut_info(i,e).xint(g,:) = p1;
        end;
      end;
      
    end;
  end;
end;

end
