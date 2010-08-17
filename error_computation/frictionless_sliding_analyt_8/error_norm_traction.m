% Calucates the L2 norm of the error in a problem with tri elements.  
% Jessica Sanders

%Initialize error
total_errorx = 0;

total_apprx = 0;


% Define 3-point Gauss quadrature (since polynom is of order 2)
gp = [-sqrt(3)/3 sqrt(3)/3];  % gauss points
gw = [1 1];               % gauss weights

% % Define 5-point Gauss quadrature
% gp = [-1/3*sqrt(5+2*sqrt(10/7)) -1/3*sqrt(5-2*sqrt(10/7)) 0 1/3*sqrt(5-2*sqrt(10/7)) 1/3*sqrt(5+2*sqrt(10/7))];  % gauss points
% gw = [(322-13*sqrt(70))/900 (322+13*sqrt(70))/900 128/225 (322+13*sqrt(70))/900 (322-13*sqrt(70))/900];               % gauss weights

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
      for i = 1:length(gp)

        % Get real coordinates of gauss points
        xn = 0.5*(1-gp(i))*p1(1)+0.5*(1+gp(i))*p2(1);
        yn = 0.5*(1-gp(i))*p1(2)+0.5*(1+gp(i))*p2(2);
    
%         % Evaluate shape functions
%         for b = 1:3     
%           % Get coorindates of area opposite node of concern
%           for m=1:3
%             jes = node(m,eleID); 
%             xes(m) = x(jes); 
%             yes(m) = y(jes);
%           end
% 
%           xes(b) = xn; yes(b) = yn;
% 
%           Area = det([[1 1 1]' xep' yep'])/2;
%           Larea = det([[1 1 1]' xes' yes'])/2;
% 
%           % Evaluate shape function for node 'b'
%           N(1,2*b-1) = N(1,2*b-1) + Larea/Area*seg_jcob*weights(g);    % First enrichment
%           N(2,2*b)   = N(2,2*b)   + Larea/Area*seg_jcob*weights(g);
%           N(1,2*b+5) = N(1,2*b+5) + Larea/Area*seg_jcob*weights(g);    % Second enrichment
%           N(2,2*b+6) = N(2,2*b+6) + Larea/Area*seg_jcob*weights(g);
%         end
% 
%         % compute derivatives of x and y wrt psi and eta
%         xdr = Ndr*xe'; ydr = Ndr*ye'; xds = Nds*xe';  yds = Nds*ye';
%         jcob = xdr*yds - xds*ydr;
% 
%         Area = det([[1 1 1]' xe' ye'])/2;

        % get normal traction at Gamma_e (depending on method for
        % enforcing the constraints at the interface
        switch IFmethod 
          case 0    % Lagrange multipliers
            lag_normal = seg_cut_info(k,e).lagmult(1);
          case 1    % penalty method
            lag_normal = seg_cut_info(k,e).ntracconv(1);
          case 2    % Nitsche's method
            lag_normal = seg_cut_info(k,e).ntracconv(1);
        end;

%         % x and y as a function of eta (for the analytical solution)
%         x_feta = N1*xe(1) + N2*xe(2) + N3*xe(3);
%         y_feta = N1*ye(1) + N2*ye(2) + N3*ye(3);

        % Get analytical solution
        anal = fless_sliding_analyt_8_traction(xn,yn);

        % Difference between analytical and numerical solutions - error!
        ex = abs(lag_normal) - abs(anal); % take absolute values, since the
                                          % sign of 'lag_normal' depends on
                                          % the enrichment, which is chosen
                                          % arbitrarily

        % square the error
        ex2 = ex^2;

        % Assemble the squared error over the element
        element_errorx = element_errorx + ex2*gw(i)*seg_jcob;

        % Assemble the approximate soln norm over the element
        approx_solnx = approx_solnx + (anal)^2*gw(i)*seg_jcob;

      end

      total_errorx = total_errorx + element_errorx;

      total_apprx = total_apprx + approx_solnx;

    end
  end
end


L2norm = sqrt(total_errorx/total_apprx);
disp(['L2-norm of normal traction:   ' num2str(L2norm)]);
