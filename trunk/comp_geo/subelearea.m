% subelearea.m
%
% Computes the area of all subelements in order to investigate some trouble
% with Nitsche's method, when single elements are cut, such that one part
% is very small.

smaller_area = [];

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
      ele_area = abs(det([[1 1 1]' xcoords' ycoords'])/2);
      
      % initialize area of part of elements that covers the positive and
      % the negative grain
      area_pos = 0;
      area_neg = 0;
      
      % loop over subelements
      for m = 1:SUBELEM_INFO(eleID).no_kids 
        % get global ID of subelement
        s = SUBELEM_INFO(eleID).kids(m);

        % get nodes of subelement
        subnodes = CONN(:,s);       

        % get coordinates of subelement nodes
        Xcoords = X(subnodes);
        Ycoords = Y(subnodes);
        
        % compute area of current subelement
        sub_area = abs(det([[1 1 1]' Xcoords' Ycoords']) / 2);
        
        if SUBELEMENT_GRAIN_MAP(s) == seg_cut_info(i,e).positive_grain
          area_pos = area_pos + sub_area;
        else
          area_neg = area_neg + sub_area;
        end;
      end;
      small_area = min(area_pos,area_neg);
      percent = small_area / ele_area;
      
      % assign the smaller area to 'smaller_area'
      smaller_area = [  smaller_area
                        seg_cut_info(i,e).elemno small_area percent];
      
    end;
  end;
end;

figure(2);
hold on;
plot(smaller_area(:,1),smaller_area(:,3),'*r');
xlabel('global element ID');
ylabel('percentage of smaller part of element based on element area');

min_percent = min(smaller_area(:,3));
disp(['Smallest fraction of area: ' num2str(min_percent)]);