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

%% first, collect data for all uncut elements
for i=1:length(x) % loop over all nodes
  % nodal coordinates
  VTKcoords = [ VTKcoords;
                x(i) y(i)];
              
  % nodal displacements
  VTKdis = [VTKdis;
            dis(2*i-1) dis(2*i)];
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
    end;
  end;
end;
% ----------------------------------------------------------------------- %

