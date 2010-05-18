% Due to the frictionless sliding, there is no tangential lagrange
% multiplier. So, only one plot will be generated.

% plot analytical solution
% figure(2);hold on;line([-1 1],[2 -2],'Color','r','LineWidth',2)

% create a new figure (no subplots due to frictionless sliding)
figure(2);      
hold on;
set(2,'Name','Lagrange multipliers (frictionless sliding)');
% title('normal direction');
ylabel('y-coordinate');
xlabel('traction value');
% axis([0 16 -1 3]);

% loop over all interfaces
for i = 1:size(seg_cut_info,1)      % every interface 'i'
  for e = 1:size(seg_cut_info,2)  % every cut element 'e' in 
                                  % interface 'i'
    if isempty(seg_cut_info(i,e).lagmult)==0    %only,if lagrange multiplier exists

      % get 2 points, that determine the subsegment
      if all(size(seg_cut_info(i,e).xint) == [2 2])       % no triple junction in 'e'
        p1 = seg_cut_info(i,e).xint(1,:);
        p2 = seg_cut_info(i,e).xint(2,:);
      elseif all(size(seg_cut_info(i,e).xint) == [1 2])   % triple junction in 'e'
        p1 = seg_cut_info(i,e).xint(1,:);

        % Get coordinates of nodes of element
        xep=zeros(1,3);
        yep=zeros(1,3);
        for m=1:3
          jep = node(m,seg_cut_info(i,e).elemno); 
          xep(m) = x(jep); 
          yep(m) = y(jep);
        end

        % Second endpoint of segment is also end point of
        % interface --> check, which one of the two endpoints
        % of the interface lies in element 'e'

        % get first endpoint
        endpoint = INTERFACE_MAP(i).endpoints(1,:); 

        % check, if it is in the element 'e'
        inside = polygon_contains_point_2d ( 3, [xep;yep], endpoint );

        if inside       % endpoint is in element
          p2 = endpoint;
        else            % endpoint is not in element
          p2 = INTERFACE_MAP(i).endpoints(2,:);
        end;
      end;
      % now, p1 and p2 are the two nodes, that determine the
      % subsegment
      xcoord = [p1(1) p2(1)];
      ycoord = [p1(2) p2(2)];

      % get lagrange multiplier for this subsegment (x and y)
      lag = seg_cut_info(i,e).lagmult;
      
      if IFmethod == 0 %|| IFmethod == 2
        % use Lagrange multipliers or Nithsche's emthod to enforce 
        % constraints at interface ('lag' is already dotted with normal)
        lag_normal = lag(1);
      else
        % use penalty or Nitsche's method
        % compute normal value
        normal = seg_cut_info(i,e).normal;
        lag_normal = lag * normal;
      end;

      % plot interface
%       line(xcoord,[lag_normal lag_normal]);
%       plot(xcoord,[lag_normal lag_normal],'-','LineWidth',3);  % horizontal interface
      plot([lag_normal lag_normal],ycoord,'-','LineWidth',3);    % vertical interface
    end;
  end;
end;  

% give a legend
% legend('analytical','simulation');