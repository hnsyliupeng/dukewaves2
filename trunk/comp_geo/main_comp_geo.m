% main_comp_geo.m
%
% mesh-generating method (based on 'polycry.m')
%

disp('generation of background mesh...');

% load data from input file 'xfeminputdata_comp_geo.mat'
load xfeminputdata_comp_geo.mat

%gets array of points p for the voronoi diagram
switch IFmeshstructure
    case 0
        structured          % call routine for structured meshing
        nonphysnodevec = [];% there are no nonphysical nodes
    case 1
        warning('MATLAB:comp_geo:meshstructure','Unstructured meshing might bring up some problems.');
        unstructured        % call routine for unstructured meshing
        nonphysnodevec = [];% there are no nonphysical nodes
    case 2
        readmeshfromGMSH    % read a mesh from a gmsh-mesh-file '*.msh'
    otherwise
        error('MATLAB:comp_geo:main_comp_geo',...
            'Unvalid mesh structure ID "IFmeshstructure" in input file');
end;

% get boundary descriptions to be able to apply NBCs via integration for
% structured meshing
if IFmeshstructure == 0 && IFneumann == 1 
  filename = fullfile('BoundaryDescriptions',IFboundarydescription);
  run(filename);
end;

% fill matrix 'p'
vdata_multi;

maxngrains = size(p,1);

%read the finite element mesh
%readmesh
%unstruct
%patch
beam_l = max(x) - min(x);
beam_h = max(y) - min(y);
gen_tol = max(beam_l,beam_h)*10e-8;
%rect
for e = 1:numele
    [nconn] = mesh_check(node(1:3,e),x,y);
    node(1:3,e) = nconn;
end

global X Y CONN;
X = x; Y = y;
CONN = node;

% x and node are reserved for the macro-elements
% X and CONN will contain info for both macroand sub elements

% set up the voronoi diagram
[vx, vy] = gen_vord(p);

% set up an interface map
% 'INTERFACE_MAP(i).endpoints' stores the x- and y-coordinates of those two
% points, which define interface 'i'.
for i=1:length(vx)
    INTERFACE_MAP(i).endpoints = [vx(:,i) vy(:,i)];
end;

%get the delaunay triangulation for the voronoi points
%this will be useful in determining which cell a point is in
tri = delaunay ( p(:,1), p(:,2) );

% plot mesh with interfaces
plotmesh

plot_vord(vx,vy,1) %plots on same figure
figure(1)
axis([0,16,-3,3])

% ----------------------------------------------------------------
% Create a new structure that has the grains associated with every
% interface

% This creates a structure interface_grains, that associates each interface
% with the two grains on either side.  

interface_grains = struct('grains',[]);

% Loop over interfaces
for i = 1:size(vx,2)
   
    distance_vec1 = zeros(1,size(p,1));
    distance_vec2 = zeros(1,size(p,1));
    xtemp1 = vx(1,i);
    ytemp1 = vy(1,i);
    
    
    
    for j = 1:size(p,1)
       
        pxtemp = p(j,1);
        pytemp = p(j,2);
        
        % Compute the distance between the two points
        
        distance_vec1(j) = sqrt((xtemp1 - pxtemp)^2 + (ytemp1 - pytemp)^2);
        
        
    end

    % Sort the distace vectors to all of the grain centers, and pick off
    % the three minima.
    [new_distance_vec, distance_vec_index] = sort(distance_vec1,'ascend');
    
    if abs(new_distance_vec(3) - new_distance_vec(2)) > gen_tol
        our_grains1 = distance_vec_index(1:2);
        interface_grains(i).grains = our_grains1;
        continue
    else
        our_grains1 = distance_vec_index(1:3);
    end
    
    
    xtemp2 = vx(2,i);
    ytemp2 = vy(2,i);
    
    for j = 1:size(p,1)
       
        pxtemp = p(j,1);
        pytemp = p(j,2);
        
        % Compute the distance between the two points
        
        distance_vec2(j) = sqrt((xtemp2 - pxtemp)^2 + (ytemp2 - pytemp)^2);
        
        
    end

    % Sort the distace vectors to all of the grain centers, and pick off
    % the three minima.
    [new_distance_vec, distance_vec_index] = sort(distance_vec2,'ascend');
    
    if abs(new_distance_vec(3) - new_distance_vec(2)) > gen_tol
        our_grains2 = distance_vec_index(1:2);
        interface_grains(i).grains = our_grains2;
        continue
    else
        our_grains2 = distance_vec_index(1:3);
    end
    
    
    
    count = 0;
    for j = 1:3
       for k = 1:3
       
           if our_grains1(j) == our_grains2(k)
           
               count = count + 1;
               temp_grains(count) = our_grains1(j);
           
           end
           
           
       
       end
    end
    
    if length(temp_grains) > 2
        error('The interface grain mapping was not unique')
    end
    
    interface_grains(i).grains = temp_grains;
    
end

% store this information into 'INTERFACE_MAP', too.
for i=1:size(interface_grains,2)
    INTERFACE_MAP(1,i).grains = interface_grains(i).grains;
end;


% clear some temporary variables
clear distance_vec1 distance_vec2 distance_vec_index temp_grains count ...
    our_grains1 our_grains2  pxtemp pytemp xtemp1 xtemp2 ytemp1 ytemp2 ...
    i new_distance_vec;

% -------------------------------------------

%loop over nodes to assign them to grains
for i=1:numnod
   k = dsearch(p(:,1), p(:,2), tri, x(i), y(i));
   nodegrainmap(i) = k;
end

% clear some temporary variables
clear k;

%figure(2)
%hold on
cutlist = zeros(1,numele);
elemgrainmap = zeros(numele,maxngrains);
%loop over elements to assign them to grains
for e=1:numele
   for i=1:3
       grainnode = nodegrainmap(node(i,e)); %find the grain associated with node i
       elemgrainmap(e,grainnode) = 1;
   end
   if (sum(elemgrainmap(e,:)) > 1)
       cutlist(e) = sum(elemgrainmap(e,:)); %add element to list of those cut by at least one grain
       %plot the element
       %plot(x(node(1:3,e)),y(node(1:3,e)))
       %plot([x(node(3,e)) x(node(1,e))],[y(node(3,e)) y(node(1,e))])
   end
end

% clear some temporary variables
clear grainnode;

%plot_vord(vx,vy,2) %plots on same figure
%axis([0,1,0,1])
%drawnow

%call partitioning routine
partition2d


% Add information to the structure seg_cut_info
for i = 1:length(vx) % Loop over all interfaces
       
    % Compute the normal to the interface
    [normal,tangent] = get_normals(vx(:,i),vy(:,i));
    
    % There are two grains associated with the interface.  
    % Define one grain as positive and one grain as negative.  
    % The normal will point outward from the positive grain.
    
    % Get first grain number 
    trial_grain = interface_grains(i).grains(1);
    
    % Get interface end points
    p1 = [vx(2,i) vy(2,i)];
    p2 = [vx(1,i) vy(1,i)];
    
    value = halfplane_contains_point_2d ( p1, p2, [p(trial_grain,1) p(trial_grain,2)] );
    
    if value
       positive_grain = trial_grain;
       negative_grain = interface_grains(i).grains(2);
    else
       positive_grain = interface_grains(i).grains(2);
       negative_grain = trial_grain;
    end
    
    for j = 1:size(seg_cut_info,2)  % Loop over cut elements

        seg_cut_info(i,j).normal = normal;
        seg_cut_info(i,j).tangent = tangent;
        seg_cut_info(i,j).grains = interface_grains(i).grains;
        seg_cut_info(i,j).positive_grain = positive_grain;
        seg_cut_info(i,j).negative_grain = negative_grain;
        seg_cut_info(i,j).interface = i;
    
    end
end

% clear some temporary variables
clear normal tangent positive_grain negative_grain grains trial_grain ...
    value p1 p2;

%assign grains to elements/subelements and plot
% figure(3)
% hold on
% plot_vord(vx,vy,3)
% axis([0,16,-3,3])

global SUBELEMENT_GRAIN_MAP;
SUBELEMENT_GRAIN_MAP = [];

%initialize
for e=1:size(CONN,2)
    SUBELEMENT_GRAIN_MAP(e) = -1;
end

for e=1:numele
   assign_grains_to_elements_recursive(p,tri,e)
end

% This is where Jessica started adding on:

numgrain = size(p,1);

global NODEINFO_ARR

%initialize
grains = [1:maxngrains;zeros(2,maxngrains)];  % Grain numbers

% 'multi_grains' is the total number of 
nodeinfo = struct('multi_grains',0,'areas',grains,'elements',[],'grain',0,...
                    'enriched','false','coords',[],'nodeno',0);


for i = 1:numnod
    nodeinfo_arr(i) = nodeinfo;
    nodeinfo_arr(i).grain = nodegrainmap(i);    % grain, in which the node resides
    nodeinfo_arr(i).coords = [x(i) y(i)];       % coordinates of the node
    nodeinfo_arr(i).nodeno = i;                 % nodenumber / node ID
end

NODEINFO_ARR = nodeinfo_arr;

% assign nodes areas associated with different grains
for e = 1:numele
    % pass the father nodes to the subroutine
    fnodes = CONN(1:3,e);   
    assign_area_to_nodes_recursive(e,fnodes);
end

% clear some temporary variables
clear fnodes;

% calulate percentages
for i = 1:numnod
    
    % count the amount of grains associated with a node
    NODEINFO_ARR(i).multi_grains = sum(NODEINFO_ARR(i).areas(2,:) ~= 0);
    if NODEINFO_ARR(i).multi_grains ~= 1
        NODEINFO_ARR(i).enriched = 'true';
    end;
    
%    if NODEINFO_ARR(i).multi_grains > 1
        for j = 1:numgrain
            NODEINFO_ARR(i).areas(3,j) = NODEINFO_ARR(i).areas(2,j)...
                /sum(NODEINFO_ARR(i).areas(2,:));
        end
%    end
end

% ----------------------------------------------------------------------- %
% THIS SECTION LOOKS FOR THE ELEMENTS, EACH NODE IS ASSIGNED TO
for i=1:numele
    ele_nodes = node(:,i)';
    for e = ele_nodes
        NODEINFO_ARR(e).elements(end+1) = i;
    end;
end;

% clear some temporary variables
clear i e ele_nodes;

% ----------------------------------------------------------------------- %

% THIS SECTION JUST DEFINES INFORMATIONAL ARRAYS ABOUT PARENT AND DAUGHTER
% ELEMENTS THAT WILL BE USED LATER

%for subelements - tells what their parent is

num_sub_elems = size(CONN,2);

global PARENTELEM_INFO; 
PARENTELEM_INFO = [];
%initialize
for e=1:num_sub_elems
    PARENTELEM_INFO(e) = -1;
end

for e=1:numele
   parentid = e;
   assign_parentelem_recursive(e,parentid)
end

% clear some temporary variables
clear parentid;

%For cut elements, list the subelements.
global SUBELEM_INFO
SUBELEM_INFO = struct('parent', 0 , 'no_kids', 0 , 'kids', []);

% initialize
for i = 1:numele
    SUBELEM_INFO(i) = struct('parent', 0 , 'no_kids', 0 , 'kids', []);
end

for i = 1:numele  % loop over parents
    if PARENTELEM_INFO(i) == -1     % if element is cut
        SUBELEM_INFO(i).parent = i;
        for j = 1:num_sub_elems   % loop over subelements
                                  % and find "daughters"
            if PARENTELEM_INFO(j) == i
                SUBELEM_INFO(i).no_kids = SUBELEM_INFO(i).no_kids + 1;
                SUBELEM_INFO(i).kids = [SUBELEM_INFO(i).kids;j];
            end   
        end
    end
end

% ----------------------------------------------------------------------- %

% INTERNAL INTERFACE SUBELEMENT PAIRS

% Here, subelements that live in the same element but different grains will
% define internal interfaces.  The data structure will be a struct with a
% parent id and pairings of subelements that define interfaces.  A pairing
% will be determined by subelements that share two nodes but are in
% different grains

% pairings will be defined by HORIZONTAL pairs of subelement numbers in the
% "pairings" arrays and the cooresponding two shared nodes in the "shared"
% array

% initialize data 

global INT_INTERFACE

INT_INTERFACE = struct('parent',0,'pairings',[],'shared',[]);

for i = 1:numele
    INT_INTERFACE(i) = struct('parent',0,'pairings',[],'shared',[]);
end

for i = 1:numele
    if cutlist(i) ~= 0
        
        INT_INTERFACE(i).parent = i;
        
        for m = 1:SUBELEM_INFO(i).no_kids  % for every subelement
            sub_elno_m = SUBELEM_INFO(i).kids(m);
            grn_m = SUBELEMENT_GRAIN_MAP(sub_elno_m); % get grain
            for n = m:SUBELEM_INFO(i).no_kids  % for every subelement
                sub_elno_n = SUBELEM_INFO(i).kids(n);
                grn_n = SUBELEMENT_GRAIN_MAP(sub_elno_n); % get grain
                flag = 0;
                if grn_m ~= grn_n
                    [flag,a,b] = shared_coord(sub_elno_m,sub_elno_n);
                end
                if flag == 1
                    INT_INTERFACE(i).pairings = [INT_INTERFACE(i).pairings;
                                                   [sub_elno_m sub_elno_n]];
                    INT_INTERFACE(i).shared = [INT_INTERFACE(i).shared;
                                                 [a b]];
                end
            end
        end
    end
end

% clear some temporary variables
clear flag grn_m grn_n sub_elno_m sub_elno_n a b;
  
% ----------------------------------------------------------------------- %
% FIND ELEMENTS THAT BUILD THE BOUNDARY OF THE DOMAIN
% 'boundary_nodes' store pairs of nodes, that build the boundary of the
% domain. Now, these node pairs have to be assigned to the elements. The
% characteristicum of a 'boundary_nodes' pair is, that it is only part of
% one single element.
% 
% If 'Ifneumann' is not defined, then set it to zero for compatibility with
% old input files
if exist('IFneumann','var') == 0,IFneumann = 0;end;

boundary_eles = [];
if IFneumann    
    for i=1:numboundele
        boundnodes = boundary_nodes(:,i);
        for j=1:numele
            elenodes = node(:,j);
            if length(intersect(boundnodes,elenodes)) == 2
                boundary_eles(i) = j;
            end;
        end;
    end;
    
    % define a structure, that stores information about the boundary
    BOUNDARY = struct('nodes',[0 0],'ele',0,'coords',[],'intsecpoint',[]);
    for i=1:numboundele
        BOUNDARY(i).nodes = boundary_nodes(:,i)';
        BOUNDARY(i).ele = boundary_eles(i);
        BOUNDARY(i).coords = [x(boundary_nodes(:,i))' y(boundary_nodes(:,i))'];
    end;
else
    % set default values, which are never used, because NBCs are given
    % nodally
    boundary_nodes = 0;
    numboundele = 0;
    BOUNDARY = 0;
end;
% ----------------------------------------------------------------------- %

% scale axes of plotted mesh
ax_x = (max(X) - min(X))/max(X);
ax_y = (max(Y) - min(Y))/max(Y);
axis([min(X)-ax_x max(X)+ax_x min(Y)-ax_y max(Y)+ax_y]);
clear ax_x ax_y;

% clear some temporary variables
clear nodeinfo_arr debug_part2d;

% SAVE DATA AS AN INPUT FILE
save my_new_mesh.mat x y node X Y CONN ELEMINFO_ARR NODEINFO_ARR...
    SUBELEMENT_GRAIN_MAP beam_h beam_l cutlist elemgrainmap maxngrains...
    nodegrainmap numele numnod p vx vy INT_INTERFACE PARENTELEM_INFO...
    SUBELEM_INFO nonphysnodevec seg_cut_info INTERFACE_MAP ...
    boundary_nodes numboundele BOUNDARY;

% print message into console
disp('Mesh generation finished. Mesh was saved to file "my_new_mesh.mat".');



