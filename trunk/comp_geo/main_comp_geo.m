% polycry.m
%
% mesh-generating method
%

% load data from input file 'xfeminputdata_comp_geo.mat'
load xfeminputdata_comp_geo.mat

%gets array of points p for the voronoi diagram
switch IFmeshstructure
    case 0
        structured
    case 1
        warning('MATLAB:comp_geo:meshstructure','Unstructured meshing might bring up some problems.');
        unstructured
    otherwise
        error('MATLAB:comp_geo:main_comp_geo',...
            'Unvalid mesh structure ID "IFmeshstructure" in input file');
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
%rect
for e = 1:numele
    [nconn] = mesh_check(node(1:3,e),x,y);
    node(1:3,e) = nconn;
end

global X Y CONN;
X = x; Y = y;
CONN = node;

%x and node are reserved for the macro-elements
%X and CONN will contain info for both macro and sub elements

%set up the voronoi diagram
[vx, vy] = gen_vord(p);

%get the delaunay triangulation for the voronoi points
%this will be useful in determining which cell a point is in
tri = delaunay ( p(:,1), p(:,2) );

plotmesh

plot_vord(vx,vy,1) %plots on same figure
figure(1)
axis([0,16,-3,3])

%loop over nodes to assign them to grains
for i=1:numnod
   k = dsearch(p(:,1), p(:,2), tri, x(i), y(i));
   nodegrainmap(i) = k;
end

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
%plot_vord(vx,vy,2) %plots on same figure
%axis([0,1,0,1])
%drawnow

%call partitioning routine
partition2d

%assign grains to elements/subelements and plot
figure(3)
hold on
plot_vord(vx,vy,3)
axis([0,16,-3,3])

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
nodeinfo = struct('multi_grains',0,'areas',grains);


for i = 1:numnod
    nodeinfo_arr(i) = nodeinfo;
end

NODEINFO_ARR = nodeinfo_arr;

% assign nodes areas associated with different grains
for e = 1:numele
    % pass the father nodes to the subroutine
    fnodes = CONN(1:3,e);   
    assign_area_to_nodes_recursive(e,fnodes);
end

% calulate percentages
for i = 1:numnod
    
    % count the amount of grains associated with a node
    NODEINFO_ARR(i).multi_grains = sum(NODEINFO_ARR(i).areas(2,:) ~= 0);
    
%    if NODEINFO_ARR(i).multi_grains > 1
        for j = 1:numgrain
            NODEINFO_ARR(i).areas(3,j) = NODEINFO_ARR(i).areas(2,j)...
                /sum(NODEINFO_ARR(i).areas(2,:));
        end
%    end
end

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
    

% ----------------------------------------------------------------------- %

% SAVE DATA AS AN INPUT FILE
save my_new_mesh.mat x y node X Y CONN ELEMINFO_ARR NODEINFO_ARR...
    SUBELEMENT_GRAIN_MAP beam_h beam_l cutlist elemgrainmap maxngrains...
    nodegrainmap numele numnod p vx vy INT_INTERFACE PARENTELEM_INFO...
    SUBELEM_INFO

% print message into console
disp('Mesh generation finished. Mesh was saved to file "my_new_mesh.mat".');



