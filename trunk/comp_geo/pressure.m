% Thick walled cylinder problem

clear

global X Y CONN
global ELEMINFO_ARR
global SUBELEMENT_GRAIN_MAP
global INT_INTERFACE
global PARENTELEM_INFO
global SUBELEM_INFO
global NODEINFO_ARR

ri = 2;
ro = 3;
rave = 2.5;

beam_l = 0;
beam_h = 1;

ndiv_th = 102;
ndiv_r  = 51;

mid_r = (1 + ndiv_r)/2;

numnod = (ndiv_r+1)*(ndiv_th+1);
numele = (ndiv_r)*(ndiv_th)*2;

cutlist = zeros(1,numele);
nodegrainmap = zeros(1,numnod);
elemgrainmap = zeros(numele,2);

SUBELEM_INFO = struct('parent',0,'no_kids',0,'kids',[]);
for i = 1:numele
    SUBELEM_INFO(i).parent = 0;
    SUBELEM_INFO(i).no_kids = 0;
    SUBELEM_INFO(i).kids = [];
end

for i = 1:(ndiv_r + 1)
    for j = 1:(ndiv_th + 1)
        r((ndiv_th+1)*(i-1)+j)  = ri + (ro-ri)*(i-1)/ndiv_r;
        th((ndiv_th+1)*(i-1)+j) = (pi/2)/(ndiv_th)*(j-1);
    end
end

for i = 1:numnod
    if r(i) > rave
        nodegrainmap(i) = 2;
    else
        nodegrainmap(i) = 1;
    end
    x(i) = r(i)*cos(th(i));
    y(i) = r(i)*sin(th(i));
end

for i = 1:ndiv_r
    for j = 1:ndiv_th
        for k = 1:2
            elem = ndiv_th*2*(i-1)+2*(j-1)+k;
            if i < mid_r
                
                elemgrainmap(elem,1) = 1;
                SUBELEMENT_GRAIN_MAP(elem) = 1;
                PARENTELEM_INFO(elem) = elem;
                
            elseif i == mid_r
                
                cutlist(elem) = 2;
                elemgrainmap(elem,1) = 1;
                elemgrainmap(elem,2) = 1;
                SUBELEMENT_GRAIN_MAP(elem) = -1;
                PARENTELEM_INFO(elem) = -1;
                
            elseif i >= mid_r
                
                elemgrainmap(elem,2) = 1;
                SUBELEMENT_GRAIN_MAP(elem) = 2;
                PARENTELEM_INFO(elem) = elem;
                
            end
            if j <= (ndiv_th/2)
                if k ==1
                    node(1,elem) = (ndiv_th+1)*(i-1) + (j-1) + 1;
                    node(2,elem) = (ndiv_th+1)*i + (j-1) + k;
                    node(3,elem) = (ndiv_th+1)*i + (j-1) + k + 1;
                elseif k == 2
                    node(1,elem) = (ndiv_th+1)*(i-1) + (j-1) + 2;
                    node(2,elem) = (ndiv_th+1)*(i-1) + (j-1) + 1;
                    node(3,elem) = (ndiv_th+1)*i + (j-1) + 2;
                end
            else
                if k ==1
                    node(1,elem) = (ndiv_th+1)*(i-1) + (j-1) + 1;
                    node(2,elem) = (ndiv_th+1)*i + (j-1) + k;
                    node(3,elem) = (ndiv_th+1)*(i-1) + (j-1) + 2;
                elseif k == 2
                    node(1,elem) = (ndiv_th+1)*(i-1) + (j-1) + 2;
                    node(2,elem) = (ndiv_th+1)*i + (j-1) + 1;
                    node(3,elem) = (ndiv_th+1)*i + (j-1) + 2;
                end
            end
        end
    end
end

p = [0 0; ro ro];

maxngrains = 2;

% Let's start the hard part:  The intersection stuff

extra_x = [];
extra_y = [];

incr = (ro - ri)/ndiv_r;

r1 = (ro + ri)/2;
r2 = (ro + ri)/2;
r3 = (ro + ri)/2 - incr/2;
r4 = (ro + ri)/2 + incr/2;
for i = 1:ndiv_th
    if i <= ndiv_th/2
        i
        th1 = (pi/2)/ndiv_th*(i-1);
        th2 = (pi/2)/ndiv_th*i;
        th3 = (pi/2)/ndiv_th*(i-1);
        th4 = (pi/2)/ndiv_th*i;
        x1 = r1*cos(th1);
        x2 = r2*cos(th2);
        x3 = r3*cos(th3);
        x4 = r4*cos(th4);
        y1 = r1*sin(th1);
        y2 = r2*sin(th2);
        y3 = r3*sin(th3);
        y4 = r4*sin(th4);
        p1 = [x1 y1];
        p2 = [x2 y2];
        q1 = [x3 y3];
        q2 = [x4 y4];
        [flg,int] = lines_exp_int_2d(p1,p2,q1,q2);
        extra_x = [extra_x x1 int(1)];
        extra_y = [extra_y y1 int(2)];
    elseif i > ndiv_th/2
        i
        th1 = (pi/2)/ndiv_th*(i-1);
        th2 = (pi/2)/ndiv_th*i;
        th3 = (pi/2)/ndiv_th*i;
        th4 = (pi/2)/ndiv_th*(i-1);
        x1 = r1*cos(th1);
        x2 = r2*cos(th2);
        x3 = r3*cos(th3);
        x4 = r4*cos(th4);
        y1 = r1*sin(th1);
        y2 = r2*sin(th2);
        y3 = r3*sin(th3);
        y4 = r4*sin(th4);
        p1 = [x1 y1];
        p2 = [x2 y2];
        q1 = [x3 y3];
        q2 = [x4 y4];
        [flg,int] = lines_exp_int_2d(p1,p2,q1,q2);
        extra_x = [extra_x x1 int(1)];
        extra_y = [extra_y y1 int(2)];
    end
end

extra_x = [extra_x 0];
extra_y = [extra_y 2.5];

X = [x extra_x];
Y = [y extra_y];

% Trianglulate subelements
num_sub_elems = 3*2*ndiv_th;

nstart = numnod/2 - ndiv_th;
estart = numele/2 - ndiv_th;

tri = zeros(3,num_sub_elems);

for i = 1:ndiv_th
    for j = 1:2
        for k = 1:3
            subel = 2*3*(i-1) + 3*(j-1) + k;
            elem = estart + 2*(i-1) + j;
            PARENTELEM_INFO(numele + subel) = elem;
            if i <= ndiv_th/2
                if (j == 1) && (k == 1)
                    SUBELEMENT_GRAIN_MAP(numele + subel) = 1;
                    tri(1,subel) = nstart + i - 1;
                    tri(2,subel) = 2*i + numnod - 1;
                    tri(3,subel) = 2*i + numnod;
                elseif (j == 1) && (k == 2)
                    SUBELEMENT_GRAIN_MAP(numele + subel) = 2;
                    tri(1,subel) = 2*i + numnod - 1;
                    tri(2,subel) = nstart + ndiv_th + i;
                    tri(3,subel) = 2*i + numnod;
                elseif (j == 1) && (k == 3)
                    SUBELEMENT_GRAIN_MAP(numele + subel) = 2;
                    tri(1,subel) = 2*i + numnod;
                    tri(2,subel) = nstart + ndiv_th + i;
                    tri(3,subel) = nstart + ndiv_th + i + 1;
                elseif (j == 2) && (k == 1)
                    SUBELEMENT_GRAIN_MAP(numele + subel) = 1;
                    tri(1,subel) = nstart + i;
                    tri(2,subel) = nstart + i - 1;
                    tri(3,subel) = 2*i + numnod;
                elseif (j == 2) && (k == 2)
                    SUBELEMENT_GRAIN_MAP(numele + subel) = 1;
                    tri(1,subel) = nstart + i;
                    tri(2,subel) = 2*i + numnod;
                    tri(3,subel) = 2*i + numnod + 1;
                elseif (j == 2) && (k == 3)
                    SUBELEMENT_GRAIN_MAP(numele + subel) = 2;
                    tri(1,subel) = 2*i + numnod + 1;
                    tri(2,subel) = 2*i + numnod;
                    tri(3,subel) = nstart + ndiv_th + i + 1;
                end
            elseif i > ndiv_th/2
                if (j == 1) && (k == 1)
                    SUBELEMENT_GRAIN_MAP(numele + subel) = 1;
                    tri(1,subel) = nstart + i;
                    tri(2,subel) = nstart + i - 1;
                    tri(3,subel) = 2*i + numnod;
                elseif (j == 1) && (k == 2)
                    SUBELEMENT_GRAIN_MAP(numele + subel) = 1;
                    tri(1,subel) = nstart + i - 1;
                    tri(2,subel) = 2*i + numnod - 1;
                    tri(3,subel) = 2*i + numnod;
                elseif (j == 1) && (k == 3)
                    SUBELEMENT_GRAIN_MAP(numele + subel) = 2;
                    tri(1,subel) = 2*i + numnod - 1;
                    tri(2,subel) = nstart + ndiv_th + i;
                    tri(3,subel) = 2*i + numnod;
                elseif (j == 2) && (k == 1)
                    SUBELEMENT_GRAIN_MAP(numele + subel) = 1;
                    tri(1,subel) = nstart + i;
                    tri(2,subel) = 2*i + numnod;
                    tri(3,subel) = 2*i + numnod + 1;
                elseif (j == 2) && (k == 2)
                    SUBELEMENT_GRAIN_MAP(numele + subel) = 2;
                    tri(1,subel) = 2*i + numnod + 1;
                    tri(2,subel) = 2*i + numnod;
                    tri(3,subel) = nstart + ndiv_th + i + 1;
                elseif (j == 2) && (k == 3)
                    SUBELEMENT_GRAIN_MAP(numele + subel) = 2;
                    tri(1,subel) = 2*i + numnod;
                    tri(2,subel) = nstart + ndiv_th + i;
                    tri(3,subel) = nstart + ndiv_th + i + 1;
                end
            end
        end
    end
end

CONN = [node tri];

num_sub_elems = size(CONN,2);

ELEMINFO_ARR = struct('nb_subelts',0,'subelemids',[]);
for i = 1:num_sub_elems
    ELEMINFO_ARR(i).nb_subelts = 1;
    ELEMINFO_ARR(i).subelemids = -1;
end

for i = 1:ndiv_th
    for j = 1:2
        elem = estart + 2*(i-1) + j;
        ELEMINFO_ARR(elem).nb_subelts = 3;
        SUBELEM_INFO(elem).parent = elem;
        SUBELEM_INFO(elem).no_kids = 3;
        for k = 1:3
            subel = 2*3*(i-1) + 3*(j-1) + k;
            subel = numele + subel;
            SUBELEM_INFO(elem).kids = [SUBELEM_INFO(elem).kids subel];
            if k == 1
                ELEMINFO_ARR(elem).subelemids = subel;
            else
                ELEMINFO_ARR(elem).subelemids = [ELEMINFO_ARR(elem).subelemids subel];
            end
        end
    end
end

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



numgrain = size(p,1);


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

% SAVE DATA AS AN INPUT FILE


 save pressure51.mat x y node X Y CONN ELEMINFO_ARR NODEINFO_ARR...
     SUBELEMENT_GRAIN_MAP cutlist elemgrainmap maxngrains...
     nodegrainmap numele numnod p INT_INTERFACE PARENTELEM_INFO...
     SUBELEM_INFO beam_h beam_l


