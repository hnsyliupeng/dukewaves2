% Two dimensional polycrystal XFEM code
% G-FEM Based on the paper by Simone, Duarte and Van der Giesseb
% All other ideas by Dolbow
% Jessica Sanders, Summer 2006 (SNL) - Summer 2007

% load input parameters from 'xfemintutdata_xfem.mat'
load xfeminputdata_xfem.mat

% assign input parameters to local variables
sliding_switch = IFsliding_switch;

% LOAD MODEL DATA
load my_new_mesh_with_BCs.mat

% ----------------------------------------------------------------------- %
    
% DETERMINE NODAL ENRICHMENTS

% determine nodes to be enriched and with which fncs
% loop over nodes, then loop over grains
% Add that grain enrichment to all of the nodes contained within

% Since each node needs to be enriched with one less grain than the total
% number which cross it's support, the algorithm is to add all of the
% enrichments possible to the data structure, and then re-loop through and
% eliminate the last enrichement.

% NODAL_ENRICH is a stucture with grain enrichments for each node, as well
% as the total number of enrichments (should be 1 or 2 in these cases).

disp('determining nodal enrichments ...');

global NODAL_ENRICH
NODAL_ENRICH = struct('cnt',0,'enrichment',[0 0 0]);

for i = 1:numnod
    NODAL_ENRICH(i) = struct('cnt',0,'enrichment',[0 0 0]);
    for j = 1:maxngrains
        % If the node belongs to a "cut" element
        if NODEINFO_ARR(i).areas(3,j) ~= 0
            
            % Keep a count of how many enrichments have been added
            NODAL_ENRICH(i).cnt = NODAL_ENRICH(i).cnt + 1; 
            cnt = NODAL_ENRICH(i).cnt;
            
            % Add the appropriate enrichment
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

% ID ARRAY FOR EQUATION NUMBERING

% Data structures for setting up the extended stiffness matrix and
% determining its size

% Here we set up an array to keep track of which nodes contain which global
% degrees of freedom.

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

disp('assembling stiffness ...');

% initialize the total number of extra added equations
multipliers = 0;

% Every cut element adds an extra set of equations

for i = 1:numele                    % for every element
    if cutlist(i) ~= 0              % if cut element
        
        multipliers = multipliers + 1;  % Count total number of multipliers
        
    end
end

old_size = numeqns;

switch sliding_switch
    case 0              % no sliding
        numeqns = numeqns + 2*multipliers;
    case 1              % frictionless sliding
        numeqns = numeqns + multipliers;
    case 2              % perfect plasticity
        warning('MATLAB:XFEM:main_xfem',...
            'There exists no code for perfect plasticity, yet.')
    case 3              % frictional contact (Coulomb)
        warning('MATLAB:XFEM:main_xfem',...
            'There exists no code for frictional contact (Coulomb), yet.')
    otherwise
        warning('MATLAB:XFEM:main_xfem',...
            'Unvalid slidingID. Choose valid ID or add additional case to switch-case-structure')
end;

% Allocate size of the stiffness matrix
%bigk = zeros(numeqns);
bigk = spalloc(numeqns,numeqns,numele*36);
ndof = 2; % number of standard degrees of freedon per node

% Get local stiffness matrices
% loop over elements

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

% APPLY CONSTRAINTS AT GRAIN INTERFACES USING LAGRANGE MULTIPLIERS

disp('enforcing constraint at interfaces ...');

ex_dofs = 0;
lag_surf = [];

for i = 1:numele  
    % for every element
    if cutlist(i) ~= 0              % if cut element
        
        % Count the extra degree of freedom we're at
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
            
            % Sliding
            
            switch sliding_switch
                case 0              % no sliding at all (fully constrained)
                case 1              % frictionless sliding
                    % Get the normal
                    [norm] = get_norm(i,j);
                case 2              % perfect plasticity
                    warning('MATLAB:XFEM:main_xfem',...
                        'There exists no code for perfect plasticity, yet.')
                case 3              % frictional contact (Coulomb)
                    warning('MATLAB:XFEM:main_xfem',...
                    'There exists no code for frictional contact (Coulomb), yet.')
                otherwise
                    warning('MATLAB:XFEM:main_xfem',...
                        'Unvalid slidingID. Choose valid ID or add additional case to switch-case-structure')
            end;
            
            % Get local constraint equations
            [ke_lag,id_node,id_lag] =...
                gen_lagrange(node,x,y,i,j,id_eqns,id_dof,pn_nodes...
                ,pos_g,neg_g,old_size,ex_dofs,sliding_switch);
            
            % If the problem includes sliding, dot with the normal
            switch sliding_switch
                case 0
                    dim = 2;
                case 1          % frictionless sliding
                    dim = 1;
                    ke_lag = ke_lag*norm';
                case 2              % perfect plasticity
                    warning('MATLAB:XFEM:main_xfem',...
                        'There exists no code for perfect plasticity, yet.')
                case 3              % frictional contact (Coulomb)
                    warning('MATLAB:XFEM:main_xfem',...
                        'There exists no code for frictional contact (Coulomb), yet.')
                otherwise
                    warning('MATLAB:XFEM:main_xfem',...
                        'Unvalid slidingID. Choose valid ID or add additional case to switch-case-structure')
            end;
            
            nlink = size(ke_lag,1);
            %
            % assemble ke_lag into bigk
            %
            for m=1:nlink
                for n=1:dim
                    rbk = id_node(m);
                    cbk = id_lag(n);
                    re = m;
                    ce = n;
                    if (rbk ~= 0) && (cbk ~= 0) 
                        % The constraint equations
                        bigk(rbk,cbk) = bigk(rbk,cbk) + ke_lag(re,ce);
                        % The transpose are the equilibrium terms
                        bigk(cbk,rbk) = bigk(cbk,rbk) + ke_lag(re,ce);
                    end
                end
            end
        end
    end
end

% ----------------------------------------------------------------------- %

% LOOP OVER "ENRICHED SURFACES" AND APPLY BCS THERE
% This is only necessary if we are applying Dirichlet boundary conditions
% on enriched nodes, and required an extensive input setup.  This works,
% but isn't well constructed.  Those boundary conditions should only be
% used if absolutely necessary.

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


% ----------------------------------------------------------------------- %
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
                    n;
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

disp('solving ...');

% Fdisp will be a vector with the solution for all global degrees of
% freedom, with "base" and enriched degrees of freedom separated.

fdisp = big_force'/bigk;

% Reassemble displacement vector  - Enriched nodes need to have their
% degrees of freedom added to end up representing a total displacement.
% The extra degrees of freedom should only be added if the node is enriched
% with the grain in which it resides.  

% ndisp(i,:) are all of the solutions (base and enriched) at node i
ndisp = zeros(numnod,6);

% dis is a vector with traditional FEM numbering, and the final solution
% for each node. It stores the nodal displacements
dis = zeros(2*numnod,1);

for i = 1:numeqns
    [nnode,doff] = find(id_eqns == i);
    ndisp(nnode,doff) = fdisp(i);
end
old_ndisp = ndisp;
for i = 1:numnod
    grain = nodegrainmap(i);
    for j = 3:6
        if id_dof(i,j) ~= grain
            ndisp(i,j) = 0;
        end
    end
    dis(2*i-1) = ndisp(i,1) + ndisp(i,3) + ndisp(i,5);
    dis(2*i) = ndisp(i,2) + ndisp(i,4) + ndisp(i,6);
end

% ----------------------------------------------------------------------- %

% POST-PROCESS

disp('postprocessing ...');

% compute stresses at center of each element
stress = zeros(numele,6,maxngrains);
strain = zeros(numele,6,maxngrains+1);
for e=1:numele
    %[stresse] = post_process(node,x,y,e,dis);
    [straine,stresse] = post_process_better(node,x,y,e,dis,old_ndisp,id_dof,cutlist,maxngrains)
    stress(e,1:6,:) = stresse;
    strain(e,1:6,:) = straine;
end



%  disp('saving to results file ...');

%  save my_results_files.mat x y node dis fdisp numele numnod cutlist...
%      INT_INTERFACE SUBELEM_INFO SUBELEMENT_GRAIN_MAP id_eqns id_dof...
%      X Y CONN stress

disp('Solving and postprocessing done.');

