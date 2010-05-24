% computemeshsize.m
%
% Computes the mesh size 'h' of structured and unstructured meshes. 'h' is defined as the
% circumradius around the biggest element.
%

% Author: Matthias Mayr (05/2010)

% initialize
h = 0;
nodes = zeros(3,1);
coords = zeros(2,3);

% computation depends on mesh structure
if IFmeshstructure == 0     % structured mesh
  % easy computation, since each elements has the same size. So, the
  % circumradius of the first element is used as the 'h'.
  
  % get global IDs of first element
  nodes = node(:,1);    
  
  % get coordinates of 'nodes' and store them in suitable format
  for i=1:3
    coords(1,i) = x(nodes(i));
    coords(2,i) = y(nodes(i));
  end;
  
  % get circumradius
  h = triangle_circumradius_2d(coords);
  
elseif IFmeshstructure == 2 %unstructured mesh (GMSH)
  % for unstructured meshes, each element can hava different size. So, you 
  % have to loop over all elements a use the biggest 'h_ele' as 'h'.
    
  % loop over all elements
  for e=1:numele
    % get global IDs of first element
    nodes = node(:,e);    

    % get coordinates of 'nodes' and store them in suitable format
    for i=1:3
      coords(1,i) = x(nodes(i));
      coords(2,i) = y(nodes(i));
    end;

    % get circumradius
    h_ele = triangle_circumradius_2d(coords);
    
    % check, if 'h_ele' is greater than 'h'
    if h < h_ele
      h = h_ele;
    end;
  end;
end;

% print mesh size 'h'
disp(['mesh size h = ' num2str(h)]);