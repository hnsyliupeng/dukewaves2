% getVTKdata.m
%
% generates the datastructures, needed to generate the VTK output
%

% Author: Matthias Mayr (08/2010)

clear VTKconn VTKcoords VTKdis VTKstress;

% initialize
VTKconn = [];
VTKcoords = [];
VTKdis = [];
VTKstress = [];
VTKcutgrid = [];

VTKinterfacenodes = [];
VTKinterfacestate = [];

%% first, collect data for all uncut elements
for i=1:length(x) % loop over all nodes
  % nodal coordinates
  VTKcoords = [ VTKcoords;
                x(i) y(i)];
              
  % nodal displacements
  VTKdis = [VTKdis;
            dis(2*i-1) dis(2*i)];
          
  VTKcutgrid = [VTKcutgrid;
                0           ];
end;

for e=1:numele        % loop over all elements 'e'
  if cutlist(e) == 0  % only, if element 'e' is not cut
    % node-element-connectivity
    VTKconn = [ VTKconn;
                node(1,e)-1 node(2,e)-1 node(3,e)-1];
    
              
    grain = SUBELEMENT_GRAIN_MAP(e);
    
    % stresses
    VTKstress = [ VTKstress;
                  stress(e,4:6,grain) stressvonmises(e,4,grain)];
  end;
end;
% ----------------------------------------------------------------------- %
%% now, collect data for all cut elements
% It might happen that some nodes
% are plotted again, but this is necessary since each node on the interface
% has to displacement values (one from each grain).
  
for e=1:numele        % loop over all elements 'e'
  if cutlist(e) ~= 0  % only, if element 'e' is cut
    % get some data from parent element 'e'
    elenodes = node(:,e);
    xcoords = x(elenodes);
    ycoords = y(elenodes);

    % loop over all subelements of element 'e'
    for m = 1:SUBELEM_INFO(e).no_kids 
      % get global ID of subelement
      s = SUBELEM_INFO(e).kids(m);

      % get grain ID of this subelement
      grain = SUBELEMENT_GRAIN_MAP(s);

      % get nodes of subelement
      subnodes = CONN(:,s);       

      % get coordinates of subelement nodes
      Xcoords = X(subnodes)';
      Ycoords = Y(subnodes)';

      curnumnodes = size(VTKcoords,1);

      % add nodes of this subelement 's'
      VTKcoords = [ VTKcoords;
                    Xcoords Ycoords];

      % add connectivity of this subelement 's'
      VTKconn = [ VTKconn;
                  curnumnodes curnumnodes+1 curnumnodes+2];

      % Interpolate displacements of current subelement's nodes
      % element displacement vector
      for b=1:3
          % base degrees of freedom
          b1 = id_eqns(elenodes(b),1);
          b2 = id_eqns(elenodes(b),2); 
          dispx(b) = totaldis(b1);
          dispy(b) = totaldis(b2);
      end

      % element displacement vector
      if id_dof(elenodes(1),3) == grain
        for b=1:3
          % extra degrees of freedom
          b1 = id_eqns(elenodes(b),3);
          b2 = id_eqns(elenodes(b),4); 
          if b1 ~= 0
            dispx(b) = dispx(b) + totaldis(b1);
            dispy(b) = dispy(b) + totaldis(b2);
          end;
        end;
      end

      if id_dof(elenodes(1),5) == grain
        for b=1:3
          % extra degrees of freedom
          b1 = id_eqns(elenodes(b),5);
          b2 = id_eqns(elenodes(b),6); 
          if b1 ~= 0
            dispx(b) = dispx(b) + totaldis(b1);
            dispy(b) = dispy(b) + totaldis(b2);
          end;
        end;
      end

      % get nodal displacements at parent nodes
      xd = dispx;
      yd = dispy;

      % compute area of parent element
      Area = det([[1 1 1]' xcoords' ycoords'])/2;

      % reset shape function matrix
      N = zeros(1,3);

      % interpolate the nodal displacements at the subelements nodes
      for n=1:3 % loop over the subelement's nodes
        xn = Xcoords(n);
        yn = Ycoords(n);

        N = zeros(1,3);

        % assemble 'N'
        for k=1:3
          % get area opposite of node of concern  
          xes = xcoords;
          yes = ycoords;
          xes(k) = xn;
          yes(k) = yn;
          Larea = det([[1 1 1]' xes' yes'])/2;

          N(k) = N(k) + Larea/Area;
        end;

        xDis(n) = N * xd';
        yDis(n) = N * yd';
      end;

      % add displacement values of subelements nodes
      VTKdis = [VTKdis;
                xDis' yDis'];    
              
      % stresses
      VTKstress = [ VTKstress;
                    stress(e,4:6,grain) stressvonmises(e,4,grain)];
                  
      % cut grid
      % find the node of the subelement, which does not coincide with a 
      % node of the parent element.
      % compute midpoints of subelemntedges
      points = zeros(3,2);
      points(1,1) = (Xcoords(1) + Xcoords(2)) / 2;  % first midpoint
      points(1,2) = (Ycoords(1) + Ycoords(2)) / 2;
      points(2,1) = (Xcoords(2) + Xcoords(3)) / 2;  % second midpoint
      points(2,2) = (Ycoords(2) + Ycoords(3)) / 2;
      points(3,1) = (Xcoords(3) + Xcoords(1)) / 2;  % third midpoint
      points(3,2) = (Ycoords(3) + Ycoords(1)) / 2;
      
      % compute distance of midpoints to parent's element edges
      distance = zeros(1,3);
      for a=1:3
        distance(a) = triangle_point_dist_2d ( [xcoords;ycoords], points(a,:) );
      end;
      
      [val id] = max(distance);
      
      gridcolor = [0;0;0];
%       if id >= 2
%         a = id-1;
%       else
%         a = 3;
%       end;
      gridcolor(id) = 1;
      
      VTKcutgrid = [VTKcutgrid;
                    gridcolor];
      
    end;
  end;
end;
% ----------------------------------------------------------------------- %
%% get data to plot the evolution of the stick slip zone
%% plot interfaces
% loop over all interfaces 'i'
for i=1:size(seg_cut_info,1)
  % loop over all cut elements 'e'
  for e=1:size(seg_cut_info,2)
    % only, if 'e' is cut by 'i'
    if seg_cut_info(i,e).elemno ~= -1
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

      % The subsegment can be divided into two parts with the same
      % length, such that each part contains one Gauss point.
      % Therefor, the coordinates of the mid point 'pm' are computed.
      pm = 0.5 * (p1 + p2);

      % Define part 1 as the connection of 'p1' and 'pm' and part 2
      % respectively. Then, the real coordinates of the gauss points
      % are needed.

      % gauss points in parameter space
      gauss = [-sqrt(3)/3 sqrt(3)/3];

      % Get real coordinates of gauss points
      xn1 = 0.5*(1-gauss(1))*p1(1)+0.5*(1+gauss(1))*p2(1);  % first gauss point
      yn1 = 0.5*(1-gauss(1))*p1(2)+0.5*(1+gauss(1))*p2(2);

      xn2 = 0.5*(1-gauss(2))*p1(1)+0.5*(1+gauss(2))*p2(1);  % second gauss point
      yn2 = 0.5*(1-gauss(2))*p1(2)+0.5*(1+gauss(2))*p2(2);

      % compute distance between 'p1' and first gauss point
      dist_1 = sqrt((p1(1) - xn1)^2 + (p1(2) - yn1)^2);

      % compute distance between 'p1' and second gauss point
      dist_2 = sqrt((p1(1) - xn2)^2 + (p1(2) - yn2)^2);

      % set color depending on current slidestate
      stylecell = {'b','r'};
      if seg_cut_info(i,e).f_trial(1) > 0
        index1 = 1;
      else
        index1 = 0;
      end;
      if seg_cut_info(i,e).f_trial(2) > 0
        index2 = 1;
      else
        index2 = 0;
      end;

      % get current normalized pseudotime of entire simulation
      timecoord = (timestep) / (length(time) - 1);

      VTKinterfacenodes = [ VTKinterfacenodes
                            p1;
                            pm;
                            pm;
                            p2];
                          
      if dist_1 < dist_2
        % first gauss point lies between 'p1' and 'pm'
        VTKinterfacestate = [ VTKinterfacestate
                              index1;
                              index2];
      else  
        % first gauss point lies between 'p2' and 'pm'
        VTKinterfacestate = [ VTKinterfacestate
                              index2;
                              index1];
      end;
    end;
  end;
end;
% ----------------------------------------------------------------------- %
