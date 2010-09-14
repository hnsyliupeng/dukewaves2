% plotgap.m
%
% Plots the tangential gap for the example from 'Simone2006'. The
% corresponding input file is 'inp_Simone2006_147_49.m'.
%

% Author: Matthias Mayr (09/2010)

disp('plottangentialgap ...');

plotdata = [];

% loop over interfaces 'i'
for i=5%1:size(seg_cut_info,1)
  % loop over cut elements 'e'
  for e=1:size(seg_cut_info,2)
    % if 'e' is cut by 'i'
    if seg_cut_info(i,e).elemno ~= -1
      %% INITIALIZE
      enrich1 = zeros(1,3);
      enrich2 = zeros(1,3);
      flg = zeros(1,6);
      % ----------------------------------------------------------------- %
      %% GET SOME ELEMENT DATA
      eleID = seg_cut_info(i,e).elemno; % element ID
      elenodes = node(:,eleID);         % node IDs
      xcoords = x(elenodes);            % x-coordinates
      ycoords = y(elenodes);            % y-coordinates
      Area = det([[1 1 1]' xcoords' ycoords'])/2; % area of element
      tangent = seg_cut_info(i,e).tangent;
      % ----------------------------------------------------------------- %
      %% ESTABLISH A SET OF FLAGS
      pos_g = seg_cut_info(i,e).positive_grain;            % positive grain
      neg_g = seg_cut_info(i,e).negative_grain;            % negative grain
      pn_nodes = get_positive_new(eleID,pos_g,neg_g); % pos. and neg. nodes
      
      % First enrichment
      for n = 1:3     % loop over nodes
        % Get the "first" enrichment of node
        enrich1(n) = id_dof(elenodes(n),3);

        if enrich1(n) == pos_g
          if (pn_nodes(n,1) == 1)
            flg(n) = 1;
          else
            flg(n) = 0;
          end
        elseif enrich1(n) == neg_g
          if (pn_nodes(n,2) == 1)
            flg(n) = -1;
          else
            flg(n) = 0;
          end        
        end
      end

      % Second Enrichment
      for n = 1:3     % loop over nodes
        % Get the "second" enrichment of nodes
        enrich2(n) = id_dof(elenodes(n),5);  

        if enrich2(n) == pos_g  % If this enrichment corresponds 
                                % to the positive grain
          if (pn_nodes(n,1) == 1)
            flg(3 + n) = 1;
          else
            flg(3 + n) = 0;
          end

        elseif enrich2(n) == neg_g
          if (pn_nodes(n,2) == 1)
            flg(3 + n) = -1;
          else
            flg(3 + n) = 0;
          end        
        end
      end
      % ----------------------------------------------------------------- %
      %% GET COORDINATES OF GAUSS POINTS
      % Gauss points in parameter space
      gauss = [-sqrt(3)/3 sqrt(3)/3];
      
      % get intersectionpoints
      intersection = seg_cut_info(i,e).xint;
      if all(size(intersection) == [2 2])
        p1 = intersection(1,:);
        p2 = intersection(2,:);
        elseif all(size(intersection) == [1 2])
        p1 = intersection(1,:);

        % Second endpoint of segment is also end point of interface
        endpoint = INTERFACE_MAP(i).endpoints(1,:);

        inside = polygon_contains_point_2d ( 3, [xcoords;ycoords], endpoint );

        if inside
          p2 = endpoint;
        else
          p2 = INTERFACE_MAP(i).endpoints(2,:);
        end
      end
      % ----------------------------------------------------------------- %
      %% LOOP OVER GAUSS POINTS
      for g = 1:length(gauss)
        % Get real coordinates of gauss points
        xn = 0.5*(1-gauss(g))*p1(1)+0.5*(1+gauss(g))*p2(1);
        yn = 0.5*(1-gauss(g))*p1(2)+0.5*(1+gauss(g))*p2(2);
      
        % evaluate gap at current gauss point
        gap = evaluate_gap_gp(xn,yn,flg,xcoords,ycoords,totaldis, ...
          id_eqns(elenodes,:),Area);
        
        gap_scalar = gap' * tangent;
        gap_pl = seg_cut_info(i,e).tgapplconv(g);
        gap_el = sign(gap_scalar) * (abs(gap_scalar) - abs(gap_pl));
        
        % vertical interface
%         figure(45);
%         hold on;
% %         plot(gap_scalar,yn,'.r');
% %         plot(gap_pl,yn,'.g');
%         plot(gap_el,yn,'.g');
        
        % get local coordinate 's'
        s = sqrt((xn-2)^2 + yn^2);
        
        % horizontal interface
%         figure(46);
%         hold on;
%         plot(s, gap_scalar,'b');
%         plot(xn,gap_pl,'.g');
%         plot(xn,gap_el,'.r');
         
        plotdata = [plotdata;
                    s gap_scalar];

      end;
      % ----------------------------------------------------------------- %
    end;
  end;
end;

plotdata_sort = sortrows(plotdata);

figure(2);
hold on;
plot(plotdata_sort(:,1),plotdata_sort(:,2),'k:','LineWidth',3);

% configure the figure
% legend('total','plastic','elastic');

disp(['method: ' num2str(IFmethod)]);

% clear some temporary variables
clear i e eleID elenodes xcoords ycoords gauss intersection p1 p2 ...
  endpoint inside xn yn gap tangent gap_scalar gap_pl gap_el pos_g ...
  neg_g pn_nodes enrich1 enrich2 flg Area;

disp('Finished.');