% Calucates the L2 norm of the error in a problem with tri elements.  
% Jessica Sanders

%Initialize error
total_errorx = 0;
total_errory = 0;

total_apprx = 0;
total_appry = 0;

% Define 2-point Gauss quadrature (since polynom is of order 2)
gp = [-sqrt(3)/3 sqrt(3)/3];  % gauss points
gw = [1 1];                   % gauss weights

% for j = 1:numele
for k=1:size(seg_cut_info,1)
  for e=1:size(seg_cut_info,2)
    if seg_cut_info(k,e).elemno ~= -1 & seg_cut_info(k,e).interface == 1

      % get global element ID
      eleID  = seg_cut_info(k,e).elemno;
      
      % get endpoints of interface
      endpoints = INTERFACE_MAP(1,seg_cut_info(k,e).interface).endpoints;

      % initialize element error
      element_errorx = 0;   % used for normal traction

      % initialize approximate solution norm
      approx_solnx = 0;

      % get nodes of parent element
      nodes = node(:,eleID);

      % get coordinates:
      xe = x(nodes);
      ye = y(nodes);
      
      % Get coordinates of parent element
      for m=1:3
        jep = node(m,eleID); 
        xep(m) = x(jep); 
        yep(m) = y(jep);
      end

      % end points of intersection - direction doesn't matter - this is for the
      % segment jacobian calculation
      intersection = seg_cut_info(k,e).xint;  % get intersection points
      if all(size(intersection) == [2 2])     % no triple junction in the element
        % get the two intersection points
        p1 = intersection(1,:);
        p2 = intersection(2,:);
      elseif all(size(intersection) == [1 2]) % triple junction in the element
        % get the two points, that define the subsegment
        p1 = intersection(1,:);

        % Second endpoint of segment is also end point of interface
        endpoint = endpoints(1,:);

        inside = polygon_contains_point_2d ( 3, [xep;yep], endpoint );

        if inside
          p2 = endpoint;
        else
          p2 = endpoints(2,:);
        end
      end
      
      % jacobian of segment to global
      he = sqrt((p1(1)-p2(1))^2 + (p1(2)-p2(2))^2);
      seg_jcob = he/2;
      
      % shape functions are not necessary here, because the discrete
      % Lagrange multipliers do not have to be multiplied with shape
      % functions, since they are constant in one subsegment of the
      % interface.      
      
      % loop over Gauss points
      for g = 1:2

        % Difference between penalty and Nitsche solution - error!
        ex = abs(seg_cut_info(k,e).ttracconv(g)) - abs(seg_cut_info_pen(k,e).ttracconv(g));

        % square the error
        ex2 = ex^2;

        % Assemble the squared error over the element
        element_errorx = element_errorx + ex2*gw(g)*seg_jcob;

        % Assemble the approximate soln norm over the element
        approx_solnx = approx_solnx + (seg_cut_info_pen(k,e).ttracconv(g))^2*gw(g)*seg_jcob;

      end

      total_errorx = total_errorx + 0.5*element_errorx;

      total_apprx = total_apprx + 0.5*approx_solnx;

    end
  end
end


L2norm = sqrt(total_errorx/total_apprx);
disp(['L2-norm of normal traction:   ' num2str(L2norm)]);
