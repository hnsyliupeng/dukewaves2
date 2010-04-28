% Two dimensional polycrystal XFEM code
% G-FEM Based on the paper by Simone, Duarte and Van der Giesseb
% All other ideas by Dolbow
% Jessica Sanders, Summer 2006 (SNL) - Summer 2007

% load input parameters from 'xfemintutdata_xfem.mat'
load xfeminputdata_xfem.mat

% Set default values to all not defined input parameters
if exist('IFmethod','var') == 0, IFmethod = 0;end;  % Lagrange multipliers
if exist('IFpenalty','var') == 0, IFpenalty = 3.0e+5;end;   % Penalty-Parameter

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

% clear some temporary variables
clear cnt i j;

% eleminate the last enrichment function
for i = 1:numnod
    last = NODAL_ENRICH(i).cnt;
    NODAL_ENRICH(i).enrichment(last) = 0;
end

% clear some temporary variables
clear last i;

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

% 'old_size' is the number of base-DOFs plus number of enriched DOFs 
old_size = numeqns;

% initialize the total number of extra added equations
multipliers = 0;

% add extra equations only for Lagrange multipliers
if IFmethod == 0
    % Every subsegment adds two extra equations
    for i = 1:size(seg_cut_info,1)     % for every interface
        for e = 1:size(seg_cut_info,2) % for every cut element in that interface
            if seg_cut_info(i,e).elemno ~= -1
                multipliers = multipliers + 1;  % Count total number of multipliers
            end
        end
    end
    
    % update number of equations 'numeqns'
    numeqns = numeqns + 2*multipliers;
end;



% add extra equations to enforce DBCs on enriched nodes via Lagrange
% multipliers
% get nodes with x- and y-DBCs
DOFs_x_DBCs = find(dispbc(1,:) ~= 0);
DOFs_y_DBCs = find(dispbc(2,:) ~= 0);

% get enriches nodes with x- and y-DBCs
enr_DOFs_x_DBCs = find(id_eqns(DOFs_x_DBCs,3) ~= 0);
enr_DOFs_y_DBCs = find(id_eqns(DOFs_y_DBCs,4) ~= 0);

enr_DOFs_x_DBCs2 = find(id_eqns(DOFs_x_DBCs,5) ~= 0);
enr_DOFs_y_DBCs2 = find(id_eqns(DOFs_y_DBCs,6) ~= 0);

% add one equation per DBC on an enriched node, but only, if this node lies
% in an enriched grain.
%
% initialize
extra_eqns_DBCx = 0;    % number of extra equations for x-DOFs as first enrichments
extra_eqns_DBCx2 = 0;   % number of extra equations for x-DOFs as second enrichments
extra_eqns_DBCy = 0;    % number of extra equations for y-DOFs as first enrichments
extra_eqns_DBCy2 = 0;   % number of extra equations for y-DOFs as second enrichments

% following variables store in following 3 values
% 1. column: global number of base DOF
% 2. column: global number of enriched DOF
% 3. column: global ID of enriched node
extra_constr_x = [];    % stores DOFs, that have to be constrained in the extra equations
extra_constr_x2 = [];   % stores DOFs, that have to be constrained in the extra equations
extra_constr_y = [];    % stores DOFs, that have to be constrained in the extra equations
extra_constr_y2 = [];   % stores DOFs, that have to be constrained in the extra equations
for i=1:length(enr_DOFs_x_DBCs)
    % get global node ID
    tempnode = DOFs_x_DBCs(enr_DOFs_x_DBCs(i));
    
    % only, if the node lies in the grain, by which it is enriched
    if nodegrainmap(tempnode) ==  id_dof(tempnode,3)
        % node lies in the grain, by which it is enriched --> add an
        % equation
        extra_eqns_DBCx = extra_eqns_DBCx + 1;
        extra_constr_x = ...
            [extra_constr_x; id_eqns(tempnode,1) id_eqns(tempnode,3) tempnode];
    end;
end;
for i=1:length(enr_DOFs_x_DBCs2)
    % get global node ID
    tempnode = DOFs_x_DBCs(enr_DOFs_x_DBCs2(i));
    
    % only, if the node lies in the grain, by which it is enriched
    if nodegrainmap(tempnode) ==  id_dof(tempnode,5)
        % node lies in the grain, by which it is enriched --> add an
        % equation
        extra_eqns_DBCx = extra_eqns_DBCx2 + 1;
        extra_constr_x2 = ...
            [extra_constr_x2; id_eqns(tempnode,1) id_eqns(tempnode,5) tempnode];
    end;
end;
for i=1:length(enr_DOFs_y_DBCs)
    % get global node ID
    tempnode = DOFs_y_DBCs(enr_DOFs_y_DBCs(i));
    
    % only, if the node lies in the grain, by which it is enriched
    if nodegrainmap(tempnode) ==  id_dof(tempnode,4)
        % node lies in the grain, by which it is enriched --> add an
        % equation
        extra_eqns_DBCy = extra_eqns_DBCy + 1;
        extra_constr_y = ...
            [extra_constr_y; id_eqns(tempnode,2) id_eqns(tempnode,4) tempnode];
    end;
end;
for i=1:length(enr_DOFs_x_DBCs2)
    % get global node ID
    tempnode = DOFs_y_DBCs(enr_DOFs_y_DBCs2(i));
    
    % only, if the node lies in the grain, by which it is enriched
    if nodegrainmap(tempnode) ==  id_dof(tempnode,6)
        % node lies in the grain, by which it is enriched --> add an
        % equation
        extra_eqns_DBCy = extra_eqns_DBCy2 + 1;
        extra_constr_y2 = ...
            [extra_constr_y2; id_eqns(tempnode,2) id_eqns(tempnode,6) tempnode];
    end;
end;

% extra_eqns_DBC = length(enr_DOFs_x_DBCs) + length(enr_DOFs_y_DBCs) + ...
%     length(enr_DOFs_x_DBCs2) + length(enr_DOFs_y_DBCs2);
extra_eqns_DBC = extra_eqns_DBCx + extra_eqns_DBCx2 + extra_eqns_DBCy ...
    + extra_eqns_DBCy2;
numeqns = numeqns + extra_eqns_DBC;



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

% clear some temporary variables
clear rbk cbk re ce nlink id ke i e j;

% ----------------------------------------------------------------------- %

% ASSEMBLE EXTERNAL FORCING VECTOR

big_force = zeros(numeqns,1);
for n = 1:numnod
    for j = 1:ndof
        slot = id_eqns(n,j);
        big_force(slot) = force(j,n);
    end
end

% clear some temporary variables
clear slot n j;

% ----------------------------------------------------------------------- %
% APPLY CONSTRAINS AT INTERFACES

% Method is chosen by paramter 'IFmethod'
switch IFmethod
    case 0          % Lagrange Multipliers
        % APPLY CONSTRAINTS AT GRAIN INTERFACES USING LAGRANGE MULTIPLIERS

        disp('enforcing constraints at interfaces via Lagrange multipliers ...');

        ex_dofs = 0;
        lag_surf = [];

         for i = 1:size(seg_cut_info,1)     % for every interface
             for e = 1:size(seg_cut_info,2) % for every cut element in that interface
                 if seg_cut_info(i,e).elemno ~= -1

                    parent_el = seg_cut_info(i,e).elemno;

        %         % Count the extra degree of freedom we're at
                    ex_dofs = ex_dofs + 1;         
                    lag_surf = [lag_surf; ex_dofs i parent_el];  % mapping between lagrange 
                                                    % multipliers and IDs of cut elements
                    ex_dofs = ex_dofs + 1;         
                    lag_surf = [lag_surf; ex_dofs i parent_el];

                    % Establish which nodes are "postively" enriched, 
                    % and which reside in the "negative" grain

                    pos_g = seg_cut_info(i,e).positive_grain;
                    neg_g = seg_cut_info(i,e).negative_grain;

                    [pn_nodes] =... 
                         get_positive_new(parent_el,pos_g,neg_g);


                     % Get local constraint equations
                     [ke_lag,id_node,id_lag] =...
                         gen_lagrange(node,x,y,parent_el,id_eqns,id_dof,pn_nodes...
                         ,pos_g,neg_g,old_size,ex_dofs,seg_cut_info(i,e).xint...
                         ,INTERFACE_MAP(i).endpoints);

                    % If the problem includes sliding, dot with the normal
                    switch sliding_switch
                        case 0
                        case 1          % frictionless sliding
                            ke_lag = ke_lag * seg_cut_info(i,e).normal;
                            ke_lag = [ke_lag zeros(12,1)];
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
                        for n=1:2
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
         end;

         % clear some temporary variables
         clear parent_el rbk cbk re ce nlink ke_lag pn_nodes i e m n id_node ...
             id_lag pos_g neg_g;

        % ---------------------------------------------------------------
        % Fix unused tangential equations
        if sliding_switch == 1
            for i = 2:2:size(lag_surf,1)
                global_index = i + max(max(id_eqns));
                bigk(global_index,global_index) = 1;
            end
        end

        % clear temporary variable 'global_index'
        clear global_index;

    case 1              % Penalty-Method
        % APPLY CONSTRAINTS AT GRAIN INTERFACES USING PENALTY METHOD
        disp('enforcing constrainst at interfaces via penalty method ...');

        % get values from input file
        penalty =IFpenalty;

        for i = 1:size(seg_cut_info,1)     % for every interface
             for e = 1:size(seg_cut_info,2) % for every cut element in that interface
                 if seg_cut_info(i,e).elemno ~= -1

                    parent_el = seg_cut_info(i,e).elemno;
                     
                    % Establish which nodes are "postively" enriched, and
                    % which reside in the "negative" grain
                    pos_g = seg_cut_info(i,e).positive_grain;
                    neg_g = seg_cut_info(i,e).negative_grain;

                    [pn_nodes] =... 
                         get_positive_new(parent_el,pos_g,neg_g);

                    % Penalty terms

                    [ke_pen,id_pen] =...
                        gen_penalty(node,x,y,parent_el,id_eqns,id_dof,...
                        pn_nodes,pos_g,neg_g,seg_cut_info(i,e).xint,...
                        INTERFACE_MAP(i).endpoints, ...
                        seg_cut_info(i,e).normal,IFsliding_switch);

                    nlink = size(id_pen,2);
                    %
                    % assemble ke_pen into bigk
                    %
                    for m=1:nlink
                        for n=1:nlink
                            rbk = id_pen(m);
                            cbk = id_pen(n);
                            re = m;
                            ce = n;
                            if ((rbk ~= 0) && (cbk ~= 0))                        
                                bigk(rbk,cbk) = bigk(rbk,cbk) + ...
                                    penalty*ke_pen(re,ce);
                            end
                        end
                    end
                end
            end
        end
        
        % clear some temporary variables
        clear rbk cbk re ce m n nlink ke_pen parent_el pos_g neg_g i e ...
            id_pen;
        
    case 2              % Nitsche's Method
        % APPLY CONSTRAINTS AT GRAIN INTERFACES USING NITSCHE'S METHOD
        disp('enforcing constraints at interfaces via Nitsche´s method ...');

        % get values from input file
        penalty =IFpenalty;

        for i = 1:size(seg_cut_info,1)     % for every interface
             for e = 1:size(seg_cut_info,2) % for every cut element in that interface
                 if seg_cut_info(i,e).elemno ~= -1

                    parent_el = seg_cut_info(i,e).elemno;
                     
                    % Establish which nodes are "postively" enriched, and
                    % which reside in the "negative" grain
                    pos_g = seg_cut_info(i,e).positive_grain;
                    neg_g = seg_cut_info(i,e).negative_grain;

                    [pn_nodes] =... 
                         get_positive_new(parent_el,pos_g,neg_g);

                    % Penalty terms

                    [ke_pen,id_pen] =...
                        gen_penalty(node,x,y,parent_el,id_eqns,id_dof,...
                        pn_nodes,pos_g,neg_g,seg_cut_info(i,e).xint,...
                        INTERFACE_MAP(i).endpoints);

                    nlink = size(id_pen,2);
                    %
                    % assemble ke_pen into bigk
                    %
                    for m=1:nlink
                        for n=1:nlink
                            rbk = id_pen(m);
                            cbk = id_pen(n);
                            re = m;
                            ce = n;
                            if ((rbk ~= 0) && (cbk ~= 0))                        
                                bigk(rbk,cbk) = bigk(rbk,cbk) + ...
                                    penalty*ke_pen(re,ce);
                            end
                        end
                    end

                    % Nitsche Terms

                    [ke_nit,id_nit] =...
                        nit_stiff(node,x,y,parent_el,id_eqns,id_dof, ...
                        pn_nodes,pos_g,neg_g,seg_cut_info(i,e).normal, ...
                        seg_cut_info(i,e).xint,INTERFACE_MAP(i).endpoints);

                    nlink = size(id_nit,2);
                    %
                    % assemble ke_nit into bigk
                    %
                    for m=1:nlink
                        for n=1:nlink
                            rbk = id_nit(m);
                            cbk = id_nit(n);
                            re = m;
                            ce = n;
                            if ((rbk ~= 0)  && (cbk ~= 0))                      
                                bigk(rbk,cbk) = bigk(rbk,cbk)...
                                    - ke_nit(re,ce) - ke_nit(ce,re);
                            end
                        end
                    end
                end
            end
        end
        
        % clear some temporary variables
        clear rbk cbk re ce m n nlink ke_pen parent_el pos_g neg_g i e ...
            id_pen ke_nit id_nit;
        
    otherwise
        error('MATLAB:XFEM:UnvalidID',...
            'Unvalid method ID. Choose valid ID or add additional case to switch-case-structure');
end;
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

% clear some temporary variables
clear cbk rbk temp re ce;


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

% clear some temporary variables
clear m m2 m3 temp;

% ----------------------------------------------------------------------- %
% ADD CONASTRAINT EQUATIONS FOR DIRICHLET BOUNDARY CONDITIONS ON ENRICHES
% NODES
% Order of assembly:
%   1. x-DBCs of first enriched nodes
%   2. y-DBCs of first enriched nodes
%   3. x-DBCs of second enriched nodes (if existing)
%   4. y-DBCs of second enriched nodes (if existing)

% x-DBCs (first enrichment)
for i=1:size(extra_constr_x,1)
    % get global base dof 'b_dof' and global enriched dof 'e_dof'
    b_dof = extra_constr_x(i,1);
    e_dof = extra_constr_x(i,2);
    
    % put in extra equation as row-equation
    rbk = old_size + 2*multipliers + i;
    bigk(rbk,b_dof) = 1.0;
    bigk(rbk,e_dof) = 1.0;
    
    % put in extra equation as column-equation (as transpose)
    bigk(b_dof,rbk) = 1.0;
    bigk(e_dof,rbk) = 1.0;
    
    % set value for prescribed displacement
    big_force(rbk) = ubar(1,extra_constr_x(i,3));
end;

% y-DBCs (first enrichment)
for i=1:size(extra_constr_y,1)
    % get global base dof 'b_dof' and global enriched dof 'e_dof'
    b_dof = extra_constr_y(i,1);
    e_dof = extra_constr_y(i,2);
    
    % put in extra equation as row-equation
    rbk = old_size + 2*multipliers + extra_eqns_DBCx + i;
    bigk(rbk,b_dof) = 1.0;
    bigk(rbk,e_dof) = 1.0;
    
    % put in extra equation as column-equation (as transpose)
    bigk(b_dof,rbk) = 1.0;
    bigk(e_dof,rbk) = 1.0;
    
    % set value for prescribed displacement
    big_force(rbk) = ubar(2,extra_constr_y(i,3));
end;

% x-DBCs (second enrichment)
for i=1:size(extra_constr_x2,1)
    % get global base dof 'b_dof' and global enriched dof 'e_dof'
    b_dof = extra_constr_x2(i,1);
    e_dof = extra_constr_x2(i,2);
    
    % put in extra equation as row-equation
    rbk = old_size + 2*multipliers + extra_eqns_DBCx + extra_eqns_DBCy + i;
    bigk(rbk,b_dof) = 1.0;
    bigk(rbk,e_dof) = 1.0;
    
    % put in extra equation as column-equation (as transpose)
    bigk(b_dof,rbk) = 1.0;
    bigk(e_dof,rbk) = 1.0;
    
    % set value for prescribed displacement
    big_force(rbk) = ubar(1,extra_constr_x2(i,3));
end;

% y-DBCs (second enrichment)
for i=1:size(extra_constr_y2,1)
    % get global base dof 'b_dof' and global enriched dof 'e_dof'
    b_dof = extra_constr_y2(i,1);
    e_dof = extra_constr_y2(i,2);
    
    % put in extra equation as row-equation
    rbk = old_size + 2*multipliers + extra_eqns_DBCx + ...
        extra_eqns_DBCx + extra_eqns_DBCy + extra_eqns_DBCx2 + i;
    bigk(rbk,b_dof) = 1.0;
    bigk(rbk,e_dof) = 1.0;
    
    % put in extra equation as column-equation (as transpose)
    bigk(b_dof,rbk) = 1.0;
    bigk(e_dof,rbk) = 1.0;
    
    % set value for prescribed displacement
    big_force(rbk) = ubar(2,extra_constr_y2(i,3));
end;

% x-DBCs (first enrichment)
% for i=1:length(enr_DOFs_x_DBCs)
%     % get global node ID
%     tempnode = DOFs_x_DBCs(enr_DOFs_x_DBCs(i));
%     
%     % only, if the node lies in the grain, by which it is enriched
%     if nodegrainmap(tempnode) ==  id_dof(tempnode,3)
%         % node lies in the grain, by which it is enriched --> add an
%         % equation
%             
%         m = id_eqns(DOFs_x_DBCs(enr_DOFs_x_DBCs(i)),3);
%         n = id_eqns(DOFs_x_DBCs(enr_DOFs_x_DBCs(i)),1);
% 
%         % base DOFs
%         bigk(old_size + 2*multipliers + i,n) = 1.0;
%         bigk(n,old_size + 2*multipliers + i) = 1.0;
% 
%         % extra DOFs
%         bigk(old_size + 2*multipliers + i,m) = 1.0;
%         bigk(m,old_size + 2*multipliers + i) = 1.0;
% 
%         % big_force
%         big_force(old_size + 2*multipliers + i) = ...
%             ubar(1,DOFs_x_DBCs(enr_DOFs_x_DBCs(i)));
%     end;
% end;

% % y-DBCs (first enrichment)
% for i=1:length(enr_DOFs_y_DBCs)
%     % get global node ID
%     tempnode = DOFs_y_DBCs(enr_DOFs_y_DBCs(i));
%     
%     % only, if the node lies in the grain, by which it is enriched
%     if nodegrainmap(tempnode) ==  id_dof(tempnode,4)
%         % node lies in the grain, by which it is enriched --> add an
%         % equation
%         m = id_eqns(DOFs_y_DBCs(enr_DOFs_y_DBCs(i)),4);
%         n = id_eqns(DOFs_y_DBCs(enr_DOFs_y_DBCs(i)),2);
% 
%         % base DOFs
%         bigk(old_size + 2*multipliers + extra_eqns_DBCx + i,n) = 1.0;
%         bigk(n,old_size + 2*multipliers + extra_eqns_DBCx + i) = 1.0;
% 
%         % extra DOFs
%         bigk(old_size + 2*multipliers + extra_eqns_DBCx + i,m) = 1.0;
%         bigk(m,old_size + 2*multipliers + extra_eqns_DBCx + i) = 1.0;
% 
%         % big_force
%         big_force(old_size + 2*multipliers + extra_eqns_DBCx + i) = ...
%             ubar(2,DOFs_y_DBCs(enr_DOFs_y_DBCs(i)));
%     end;
% end;
% 
% % x-DBCs (second enrichment)
% for i=1:length(enr_DOFs_x_DBCs2)
%     % get global node ID
%     tempnode = DOFs_x_DBCs(enr_DOFs_x_DBCs2(i));
%     
%     % only, if the node lies in the grain, by which it is enriched
%     if nodegrainmap(tempnode) ==  id_dof(tempnode,5)
%         % node lies in the grain, by which it is enriched --> add an
%         % equation
%             
%         m = id_eqns(DOFs_x_DBCs(enr_DOFs_x_DBCs2(i)),5);
%         n = id_eqns(DOFs_x_DBCs(enr_DOFs_x_DBCs2(i)),1);
% 
%         % base DOFs
%         bigk(old_size + 2*multipliers + i,n) = 1.0;
%         bigk(n,old_size + 2*multipliers + i) = 1.0;
% 
%         % extra DOFs
%         bigk(old_size + 2*multipliers + extra_eqns_DBCx + ...
%             extra_eqns_DBCy + i,m) = 1.0;
%         bigk(m,old_size + 2*multipliers + extra_eqns_DBCx + ...
%             extra_eqns_DBCy+ i) = 1.0;
% 
%         % big_force
%         big_force(old_size + 2*multipliers + extra_eqns_DBCx + ...
%             extra_eqns_DBCy + i) = ubar(1,DOFs_x_DBCs(enr_DOFs_x_DBCs(i)));
%     end;
% end;
% 
% % y-DBCs (second enrichment)
% for i=1:length(enr_DOFs_y_DBCs2)
%     % get global node ID
%     tempnode = DOFs_y_DBCs(enr_DOFs_y_DBCs2(i));
%     
%     % only, if the node lies in the grain, by which it is enriched
%     if nodegrainmap(tempnode) ==  id_dof(tempnode,6)
%         % node lies in the grain, by which it is enriched --> add an
%         % equation
%         m = id_eqns(DOFs_y_DBCs(enr_DOFs_y_DBCs2(i)),6);
%         n = id_eqns(DOFs_y_DBCs(enr_DOFs_y_DBCs2(i)),2);
% 
%         % base DOFs
%         bigk(old_size + 2*multipliers + extra_eqns_DBCx + ...
%             extra_eqns_DBCy + extra_eqns_DBCx2 + i,n) = 1.0;
%         bigk(n,old_size + 2*multipliers + extra_eqns_DBCx + ...
%             extra_eqns_DBCy + extra_eqns_DBCx2 + i) = 1.0;
% 
%         % extra DOFs
%         bigk(old_size + 2*multipliers + extra_eqns_DBCx + ...
%             extra_eqns_DBCy + extra_eqns_DBCx2 + i,m) = 1.0;
%         bigk(m,old_size + 2*multipliers + extra_eqns_DBCx + ...
%             extra_eqns_DBCy + extra_eqns_DBCx2+ i) = 1.0;
% 
%         % big_force
%         big_force(old_size + 2*multipliers + extra_eqns_DBCx + ...
%             extra_eqns_DBCy + extra_eqns_DBCx2 + i) = ...
%             ubar(2,DOFs_y_DBCs(enr_DOFs_y_DBCs(i)));
%     end;
% end;


% ----------------------------------------------------------------------- %
% FIX NONPHYSICAL NODES (due to gmsh-meshes)

for i = nonphysnodevec
    dofvec_temp = id_eqns(i,:); % get global DOF-numbers for nonphysical node
    for j=dofvec_temp
        bigk(j,j)=1;
    end;
end;

% clear some temporary variables
clear dofvec_temp;


% ----------------------------------------------------------------------- %

% LINEAR SOLVE AND RE-ASSEMBLE SOLUTION 

% solve stiffness equations

disp('solving ...');

% Fdisp will be a vector with the solution for all global degrees of
% freedom, with "base" and enriched degrees of freedom separated.

switch IFSolverType
    case 0                          % explicit solver
        disp('    Solver: explicit');
        disp(['    Condition number in L1-norm:  ' num2str(condest(bigk))]);
        fdisp = big_force'/bigk;    % fdisp is a row-vector
    case 1                          % implicit solver (Newton-scheme)
        disp('Solver: implicit');
        
        % implizit solving via a Newton-Raphson-Scheme
        % set some default values (if they are not set via input file)
        if IFmaxiter==1;maxiter = 25;end;  % set 'maxiter=25' (default)
        if IFconvtol==0;convtol = 1.0e-5;end;% set 'convtols=1e-5' (default)
        
        % initialize some variables
        solu = zeros(size(big_force));            % set a start vector
        iter = 0;                                 % iteration index for newton-scheme
        
        % store old solution to enable update
        old_solu = solu;
        
        while 1         % maximum number of iterations = maxiter
            % add '1' to the iteration index
            iter = iter + 1;

            % compute the increment vector 'delta'
            delta = -old_solu + bigk\big_force; 

            % update the solution-vector
            solu = solu + delta;

            % print some information about the current iteration step
            step_info = ['  Newton step: ' num2str(iter) '    Res-Norm: '...
                num2str(norm(delta))];
            disp(step_info);
            
            % convergence check
            if norm(delta) < IFconvtol;break;end;

            % save displacement to 'old_solu'
            old_solu = solu;

            if iter > IFmaxiter
                error('MATLAB:XFEM:main_xfem',...
                    'Newton did not converge in %d iterations.', IFmaxiter);
            end;
        end;
        fdisp = solu';      % fdisp is a row-vector
        
        % clear some temporary variables
        clear solu step_info delta old_solu iter;
    otherwise
        error('MATLAB:XFEM:main_xfem',...
            'Unvalid solver type ID. Choose a valid ID or add an additional solver in "main_xfem.m".');
end;

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

% computing the internal forces in the interface depends on the method of
% constraint enforcing
switch IFmethod 
    case 0              % Lagrange multipliers
        % extract vector with lagrange multipliers from 'fdisp'
        lagmult = fdisp(end-2*multipliers+1:end);

        % assign lagrange multipliers into 'seg_cut_info'
        % a subsegment is defined uniquely by interface and element ID
        % loop over interfaces 'i' and elements 'e'

        for i = 1:size(seg_cut_info,1)      % every interface 'i'
            for e = 1:size(seg_cut_info,2)  % every element 'e'
                % initialize variable for lagrange multiplier in 'seg_cut_info'
                seg_cut_info(i,e).lagmult = [];
                
                if isempty(seg_cut_info(i,e).elemno)==0   % only, if element is cut by 'i'
                    eleID = seg_cut_info(i,e).elemno;

                    % find row in 'lagmult', that suits to current interface 'e'
                    % and to the element, that is cut by 'i'
                    index = find(lag_surf(:,2)==i & lag_surf(:,3)==eleID);
                    seg_cut_info(i,e).lagmult = lagmult(index);
                 end;
            end;
        end;

        % clear some temporary variables
        clear eleID index i e;
        
    case 1                  % Penalty method
        % Compute Lagrange multipliers via 'lambda = alpha * [[u]]'
        
        %initialize
        lagmult = [];
        
        for i = 1:size(seg_cut_info,1)     % for every interface
            for e = 1:size(seg_cut_info,2) % for every cut element in that interface
                % initialize variable for lagrange multiplier in 'seg_cut_info'
                seg_cut_info(i,e).lagmult = [];
                
                if seg_cut_info(i,e).elemno ~= -1

                    parent_el = seg_cut_info(i,e).elemno;

                    % Establish which nodes are "postively" enriched, and
                    % which reside in the "negative" grain
                    pos_g = seg_cut_info(i,e).positive_grain;
                    neg_g = seg_cut_info(i,e).negative_grain;

                    [pn_nodes] =... 
                         get_positive_new(parent_el,pos_g,neg_g);
        
                     seg_cut_info(i,e).lagmult = ...
                         get_lag_mults_for_penalty(node,x,y,parent_el, ...
                         id_eqns,id_dof,pn_nodes,pos_g,neg_g, ...
                         seg_cut_info(i,e).xint, ...
                         INTERFACE_MAP(i).endpoints,penalty,fdisp);
                     
                     lagmult = [lagmult seg_cut_info(i,e).lagmult];
                end;
            end;
        end;
                     
        % clear some temporary variables
        clear pn_nodes neg_g pos_g i e parent_el;
          
    case 2                  % Nitsche's method
    otherwise
        error('MATLAB:XFEM:UnvalidID',...
            'Unvalid method ID. Choose a valid ID or add an additional case to switch-case-structure.');
end;
% ----------------------------------------------------------------------- %

% POST-PROCESS

disp('postprocessing ...');

% compute stresses 'stress' and strains 'strain' at center of each element
% Structure of 'stress':
%   Dimension: numelex6x3
%   Index 1:    'stress' for element 'e'
%   Index 2:    column 1    global element ID
%               column 2    x-coordinate of element centroid
%               column 3    y-coordinate of element centroid
%               column 4    xx-stress at element centroid
%               column 5    yy-stress at element centroid
%               column 6    xy-stress at element centroid
%   Index 3:    ID of grain, to which these values belong to
stress = zeros(numele,6,maxngrains);
% 'strain' has an equivalent structure
strain = zeros(numele,6,maxngrains+1);
for e=1:numele
%     [stresse] = post_process(node,x,y,e,dis);
%     stress(e,1:6) = stresse;
        
    [straine,stresse] = post_process_better(node,x,y,e,dis,old_ndisp,id_dof,cutlist,maxngrains);
    stress(e,1:6,:) = stresse;
    strain(e,1:6,:) = straine;
end

% clear some temporary variables
clear stresse straine maxstress_vec minstress_vec i e f j;


%  disp('saving to results file ...');

%  save my_results_files.mat x y node dis fdisp numele numnod cutlist...
%      INT_INTERFACE SUBELEM_INFO SUBELEMENT_GRAIN_MAP id_eqns id_dof...
%      X Y CONN stress

disp('Solving and postprocessing done.');

