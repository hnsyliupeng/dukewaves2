% Two dimensional polycrystal XFEM code
% G-FEM Based on the paper by Simone, Duarte and Van der Giesseb
% All other ideas by Dolbow
% Jessica Sanders, Summer 2006 (SNL) - Summer 2007

%clear

% ----------------------------------------------------------------------- %

% LOAD MODEL DATA

sprintf('loading model data')

load cstruct4.mat
%load unstruct2.mat
% load pressure51.mat

num_sub_elems = size(CONN,2);

% ----------------------------------------------------------------------- %

% GRAIN MATERIAL INFORMATION

global GRAININFO_ARR
GRAININFO_ARR = struct('grain_no', 0, 'num_elems',0,'poisson',0,'youngs',0);

% assign material properties to each grain
poissons = [0.3 0.3 0.3 0.3];
youngs = [1000.0 1000.0 1000.0 1000.0];
%youngs = [3*10^7 3*10^7 3*10^7 3*10^7];
neg = sum(elemgrainmap);    % number of elements in each grain

for i = 1:maxngrains
    GRAININFO_ARR(i).grain_no = i;
    GRAININFO_ARR(i).num_elems = neg(i);
    GRAININFO_ARR(i).poisson = poissons(i);
    GRAININFO_ARR(i).youngs = youngs(i);
end

% ----------------------------------------------------------------------- %
    
% DETERMINE NODAL ENRICHMENTS

% determine nodes to be enriched and with which fncs
% loop over grains, adding that grain enrichment to all of the elements
% contained within

sprintf('determining nodal enrichments')

global NODAL_ENRICH
NODAL_ENRICH = struct('cnt',0,'enrichment',[0 0 0]);

for i = 1:numnod
    NODAL_ENRICH(i) = struct('cnt',0,'enrichment',[0 0 0]);
%    for j = maxngrains:-1:1
    for j = 1:maxngrains
        if NODEINFO_ARR(i).areas(3,j) ~= 0
            NODAL_ENRICH(i).cnt = NODAL_ENRICH(i).cnt + 1; 
            cnt = NODAL_ENRICH(i).cnt;
            NODAL_ENRICH(i).enrichment(cnt) = j;         
        end
    end
end

% eleminate the last enrichment function
for i = 1:numnod
    last = NODAL_ENRICH(i).cnt;
    NODAL_ENRICH(i).enrichment(last) = 0;
end 

% ----------------------------------------------------------------------- %
   
% BOUNDARY CONDITIONS  - FORCE AND DISPLACEMENTS

sprintf('enforcing boundary conditions')

f = 0;
[force, dispbc, ubar, num_enr_surf, enr_surf, bc_enr]...
    = applybcs(x,y,numnod,beam_l,beam_h,f);

% ----------------------------------------------------------------------- %

% Data structures for setting up the extended stiffness matrix and
% determining its size

% ID array including extra degrees of freedom
% preallocate 6 slots for each node - x = -1,y = -2,and 2 enrichments with
% two degrees of freedom each (this allows up to three way intersectons). 

id_dof = zeros(numnod,6);
id_eqns = zeros(numnod,6);

eqn_num = 1;
for i = 1:numnod
    id_dof(i,1) = -1;          % indicates x degree of freedom
    id_eqns(i,1) = eqn_num;
    eqn_num = eqn_num + 1;
    id_dof(i,2) = -2;          % indicates y degree of freedom
    id_eqns(i,2) = eqn_num;
    eqn_num = eqn_num + 1; 
    if NODAL_ENRICH(i).enrichment(1) ~= 0
        id_dof(i,3) = NODAL_ENRICH(i).enrichment(1);
        id_eqns(i,3) = eqn_num;
        eqn_num = eqn_num + 1;
        id_dof(i,4) = NODAL_ENRICH(i).enrichment(1);
        id_eqns(i,4) = eqn_num;
        eqn_num = eqn_num + 1;
    end
    if NODAL_ENRICH(i).enrichment(2) ~= 0
        id_dof(i,5) = NODAL_ENRICH(i).enrichment(2);
        id_eqns(i,5) = eqn_num;
        eqn_num = eqn_num + 1;
        id_dof(i,6) = NODAL_ENRICH(i).enrichment(2);
        id_eqns(i,6) = eqn_num;
        eqn_num = eqn_num + 1;
    end
end

numeqns = max(max(id_eqns));

% ----------------------------------------------------------------------- %


%ASSEMBLY OF STIFFNESS

sprintf('assembling stiffness')

multipliers = 0;

for i = 1:numele                    % for every element
    if cutlist(i) ~= 0              % if cut element
        
        multipliers = multipliers + 1;  % Count total number of multipliers
        
    end
end

old_size = numeqns;

numeqns = numeqns + 2*multipliers;

%bigk = zeros(numeqns);
bigk = spalloc(numeqns,numeqns,numele*36);
ndof = 2; % number of standard degrees of freedon per node
%
% loop over elements
%
for e = 1:numele
    if (ELEMINFO_ARR(e).nb_subelts == 1)
       [id,ke] = elemstiff_class(node,x,y,e,id_dof,id_eqns);
    else %call the recursive routine
       [id,ke] = elemstiff_class_recursive(node,x,y,e,e,id_dof,id_eqns);
    end
    nlink = size(id,2);

   %
   % assemble ke into bigk
   %
   for i=1:nlink;
      for j=1:nlink;
         rbk = id(i);
         cbk = id(j);
         re = i;
         ce = j;
         bigk(rbk,cbk) = bigk(rbk,cbk) + ke(re,ce);
      end
   end
end

% ----------------------------------------------------------------------- %

% ASSEMBLE EXTERNAL FORCING VECTOR

big_force = zeros(numeqns,1);
for n = 1:numnod
    for j = 1:ndof
        slot = id_eqns(n,j);
        big_force(slot) = force(j,n);
    end
end

% ----------------------------------------------------------------------- %

% APPLY CONSTRAINTS AT GRAIN INTERFACES

sprintf('enforcing constraint at interfaces')

ex_dofs = 0;
lag_surf = [];
tol = 0.00001;

for i = 1:numele  
    % for every element
    if cutlist(i) ~= 0              % if cut element
        
        ex_dofs = ex_dofs + 1;
        
        lag_surf = [lag_surf; ex_dofs i];
        
        for j = 1:size(INT_INTERFACE(i).pairings,1) % for every subsegment
                                                    % pairing
                                                    % i
                                                    % j

            % Establish which grain is "positive" and which is "negative"
            % as well as which nodes are "positively" and which are
            % "negatively" enriched. 
            
            [pos_g,neg_g,pn_nodes] =... 
                get_positive(i,j,nodegrainmap,p);
            
            
            [ke_lag,id_node,id_lag] =...
                gen_lagrange(node,x,y,i,j,id_eqns,id_dof,pn_nodes...
                ,pos_g,neg_g,old_size,ex_dofs)
        
            nlink = size(ke_lag,1);
            %
            % assemble ke_pen into bigk
            %
            for m=1:nlink
                for n=1:2
                    rbk = id_node(m);
                    cbk = id_lag(n);
                    re = m;
                    ce = n;
                    if (rbk ~= 0) && (cbk ~= 0) 
                    bigk(rbk,cbk) = bigk(rbk,cbk) + ke_lag(re,ce);
                    bigk(cbk,rbk) = bigk(cbk,rbk) + ke_lag(re,ce);
                    end
                end
            end
        end
    end
end

% ----------------------------------------------------------------------- %

% LOOP OVER "ENRICHED SURFACES" AND APPLY BCS THERE

if (bc_enr)         % The switch indicating we have enriched nodes with bcs
    
    temp = size(bigk,2);

    bigk(temp+2*num_enr_surf,1) = 0;            % Resize bigk
    bigk(1,temp+2*num_enr_surf) = 0;

    for i = 1:num_enr_surf                      % Loop over sub surfaces
        
        
        [ke_con,id_node,id_lag,fe_con] =...
        enr_constraints(node,x,y,i,id_eqns,id_dof,enr_surf,temp,ubar,dispbc,nodegrainmap);

        ke_con = ke_con';

        nlink = size(ke_con,1);
        %
        % assemble ke_pen into bigk
        %
        for m=1:nlink
            for n=1:2
                rbk = id_node(m);
                cbk = id_lag(n);
                re = m;
                ce = n;
                if (rbk ~= 0) && (cbk ~= 0) 
                    bigk(rbk,cbk) = bigk(rbk,cbk) + ke_con(re,ce);
                    bigk(cbk,rbk) = bigk(cbk,rbk) + ke_con(re,ce);
                end
            end
        end
        %
        % assemble force into RHS
        %
        big_force(id_lag(1)) = fe_con(1);
        big_force(id_lag(2)) = fe_con(2);         
    end
    for i = num_enr_surf:-1:1    
        if enr_surf(i).xy(2) ~= 1
            bigk(temp+2*i,:) = [];
            bigk(:,temp+2*i) = [];
            big_force(temp + 2*i) = [];
        end
        if enr_surf(i).xy(1) ~= 1
            bigk(temp+2*i-1,:) = [];
            bigk(:,temp+2*i-1) = [];
            big_force(temp + 2*i-1) = [];
        end
    end
end

% ENFORCE DISPLACEMENT BOUNDARY CONDITIONS 

for n=1:numnod
    for j=1:ndof
        if (dispbc(j,n) == 1)          
            m  = id_eqns(n,j);
            m2 = id_eqns(n,j+2);
            m3 = id_eqns(n,j+4);
            if (m2 == 0) && (m3 == 0)       % If the node is unenriched
                temp = size(bigk,2);
                bigk(m,:) = zeros(1,temp);
                big_force = big_force - (ubar(j,n)*bigk(:,m));
                bigk(:,m) = zeros(temp,1);
                bigk(m,m) = 1.0;
                big_force(m) = ubar(j,n);
            elseif (m2 ~= 0) && (m3 == 0)               % If the node is enriched
                if id_dof(n,j+2) ~= nodegrainmap(n);    % But only 1 dof is active
                    n
                    temp = size(bigk,2);
                    bigk(m,:) = zeros(1,temp);
                    big_force = big_force - (ubar(j,n)*bigk(:,m));
                    bigk(:,m) = zeros(temp,1);
                    bigk(m,m) = 1.0;
                    big_force(m) = ubar(j,n);
                end
            end
        end
    end
end
  
% ----------------------------------------------------------------------- %

% LINEAR SOLVE AND RE-ASSEMBLE SOLUTION 

% solve stiffness equations

sprintf('solving')

fdisp = big_force'/bigk;

% Reassemble displacement vector

ndisp = zeros(numnod,6);
disp = zeros(2*numnod,1);
for i = 1:numeqns
    [nnode,doff] = find(id_eqns == i);
    ndisp(nnode,doff) = fdisp(i);
end
for i = 1:numnod
    grain = nodegrainmap(i);
    for j = 3:6
        if id_dof(i,j) ~= grain
            ndisp(i,j) = 0;
        end
    end
    disp(2*i-1) = ndisp(i,1) + ndisp(i,3) + ndisp(i,5);
    disp(2*i) = ndisp(i,2) + ndisp(i,4) + ndisp(i,6);
end

% ----------------------------------------------------------------------- %

% POST-PROCESS

sprintf('postprocessing')

% compute stresses at center of each element

stress = zeros(numele,6);
for e=1:numele
    [stresse] = post_process(node,x,y,e,disp);
    stress(e,1:6) = stresse;
end

sprintf('saving to results file')

%  save nlag_pres51.mat x y node disp fdisp numele numnod cutlist...
%      INT_INTERFACE SUBELEM_INFO SUBELEMENT_GRAIN_MAP id_eqns id_dof...
%      X Y CONN


