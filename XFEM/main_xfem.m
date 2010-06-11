%% Two dimensional polycrystal XFEM code
% G-FEM Based on the paper by Simone, Duarte and Van der Giesseb
% All other ideas by Dolbow
% Jessica Sanders, Summer 2006 (SNL) - Summer 2007
%% LOAD INPUT PARAMETERS FROM 'xfemintutdata_xfem.mat'
load xfeminputdata_xfem.mat
% ----------------------------------------------------------------------- %
%% LOAD MODEL DATA
load my_new_mesh_with_BCs.mat
% ----------------------------------------------------------------------- %
%% SET DAFAULT VALUES TO ALL NOT DEFINED INPUT PARAMETERS
% Some variables, which are necessary for all sliding cases

% parameters for enforcing the constraints at the interface
if exist('IFmethod','var') == 0, IFmethod = 0;end;          % Lagrange multipliers
if exist('IFpenalty','var') == 0, IFpenalty = 3.0e+5;end;   % Penalty-Parameter
if exist('IFnitsche','var') == 0, IFnitsche = 3.0e+5;end;   % Stabilization-Parameter

% parameters for the loadstepping scheme
if exist('IFtime','var') == 0, IFtime = [0 1];end;          % load steps

% parameters for Quasi-Newton-scheme
if IFmaxiter==1;maxiter = 25;end;     % maximum number of iterations
if IFconvtol==0;convtol = 1.0e-5;end; % convergence tolerance

% some variables, that are necessary for specific sliding cases only
switch IFsliding_switch
  case 0  % fully tied case
  case 1  % frictionless sliding
  case 2  % perfect plasticity
    if exist('IFyieldstress','var') == 0, IFyieldstress = 1.0e+12; end; % yield stress
  case 3  % frictional sliding (Coulomb)
end;
% ----------------------------------------------------------------------- %
%% ASSIGN INPUT PARAMETERS TO LOCAL VARIABLES
% sliding parameters
sliding_switch = IFsliding_switch;

% time / load stepping parameters & displacement increments for each
% loadstep
if length(IFtime) > 1
  time = IFtime;
  dis_increment = zeros(2*numnod,1,(length(time)-1));
else
  time = [0 1];
  dis_increment = zeros(2*numnod,1,1);
end;
% ----------------------------------------------------------------------- %
%% INITIALIZE (1)
% flag for changes in 'slidestate'-flags
slidestateconv = [];

% diplacement vector with separated base and enriched DOFs before
% re-assembling
old_ndisp = zeros(numnod,6);

% vector for current displacements during the solution process
current_dis = zeros(1,2*numnod);

% counter 
count = 0;

% vector for 'Lagrange multipliers' (= tractions at interface)
lagmult = [];

% initialize sliding state flags 'slidestate' in 'seg_cut_info'
%   flag    description
%   0       stick
%   1       slip
for i=1:size(seg_cut_info,1)
  for e=1:size(seg_cut_info,2)
    seg_cut_info(i,e).slidestate = 0; % initialize as 'stick'
  end;
end;
% ----------------------------------------------------------------------- %
%% DETERMINE NODAL ENRICHMENTS

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
%% ID ARRAY FOR EQUATION NUMBERING

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
%% PREPARE ASSEMBLY OF STIFFNESS

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

% ----------------------------------------------------------------------- %
%% ASSEMBLE ELASTIC CONTRIBUTION TO GOBAL STIFFNESS MATRIX 'bigk'
% Get local stiffness matrices
% loop over elements

for e = 1:numele
  if (ELEMINFO_ARR(e).nb_subelts == 1)
   [id,ke] = elemstiff_class(node,x,y,e,id_dof,id_eqns);
  else  % call the recursive routine
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

% store elastic contribution to stiffnes matrix for later use
bigk_el = bigk;

% clear some temporary variables
clear rbk cbk re ce nlink id ke i e j;

% ----------------------------------------------------------------------- %
%% PREPARE ASSEMBLY OF EXTERNAL FORCE VECTOR
% % get enriched nodes with NBCs on it
% DOFs_x_NBCs = find(force(1,:));     % nodes with NBC in x-direction
% DOFs_y_NBCs = find(force(2,:));     % nodes with NBC in y-direction
% 
% % get enriched dofs of those nodes
% % get enriches nodes with x- and y-DBCs
% % index of those nodes in DOFs_x_NBCs, that have a first enrichment
% enr_DOFs_x_NBCs = find(id_eqns(DOFs_x_NBCs,3) ~= 0); 
% % index of those nodes in DOFs_y_NBCs, that have a first enrichment
% enr_DOFs_y_NBCs = find(id_eqns(DOFs_y_NBCs,4) ~= 0);    
% % index of those nodes in DOFs_y_NBCs, that have a second enrichment
% enr_DOFs_x_NBCs2 = find(id_eqns(DOFs_x_NBCs,5) ~= 0);   
% % index of those nodes in DOFs_y_NBCs, that have a second enrichment
% enr_DOFs_y_NBCs2 = find(id_eqns(DOFs_y_NBCs,6) ~= 0);   

% For not enriched nodes, the nodal force given in 'force' can be assigned
% to 'big_force' at the right position. For enriched nodes, the element
% load vector has to be computed via integration, because the nodal forces
% have to be partitioned on the base and extra DOFs.
%
% There are two methods to impose NBCs:
% (1)   Nodal forces, given in input file (works only on not-enriched nodes
% (2)   Give traction as a function and get nodal forces via integration
%       over the Neumann boundary.
%
% Set standard method to "nodal forces", if not specified in input file
if exist('IFneumann','var') == 0,IFneumann = 0; end;

% initialize global force vector
big_force = zeros(numeqns,1);

% initialize global vector for tractions at the interface (used to assemble
% force contributions due to tractions, caused by plasticity or frictional 
% sliding at the interface. 
big_force_traction = zeros(numeqns,1);
% ----------------------------------------------------------------------- %
%% ASSEMBLE EXTERNAL GLOBAL FORCE VECTOR 'big_force'
% select method of applying Neumann BCs
switch IFneumann
  case 0  % give nodal forces in input file
    for n = 1:numnod
      for j = 1:ndof
        slot = id_eqns(n,j);
        big_force(slot) = force(j,n);
      end

            % If node is enriched, then apply NBCs on enriched DOFs, too.
%             if strcmp(NODEINFO_ARR(n).enriched,'true')
%                 % These nodes are enriched, so the forces can not be applied
%                 % directly to the DOFs. The integration over Gamma_h has to be
%                 % split into 2 parts at the interface.
%                 %
%                 % Input File: Give triples of nodes, that contain the enriched node
%                 % with NBC (in first column) and its two neighbours in ascending 
%                 % order (in columns 2 and 3), in variable 'nodeNBC'
% 
%                 % apply NBCs only on nodes, that are element of the Neumann
%                 % boundary. For eniched nodes, these nodes have to be listed in
%                 % 'nodeNBC' in the NBC-input-file by the user (together with the
%                 % two neighbored nodes).
%                 if any(n == nodeNBC(:,1))
%         %             disp(['enriched node n = ' num2str(n)]);
%                     % get the two nodes (in ascending order), neighbored to 
%                     % enriched node 'n'
%         %             othernodes = nodeNBC(find(nodeNBC(:,1)==n),2:3);
% 
%                     % call fucntion, that manages the integration over the Neumann
%                     % boundary
%         %             [force_values, force_id] = ...
%         %                NBCs_on_enr_nodes_new(NODEINFO_ARR([n othernodes]), ...
%         %                force(:,[n othernodes]), seg_cut_info,interface_grains, ...
%         %                id_dof,node);
% 
%                 end;
%             end;
    end;

  case 1  % give function for forces in input file
    % loop over all boundary elements
    for i=1:size(BOUNDARY,2);
    % loop over all traction-sets
      for j=1:size(FORCE,2)
        % integrate only over elements, where both nodes have
        % forces on it
        if length(intersect(BOUNDARY(1,i).nodes,FORCE(1,j).nodes)) == 2
          % check, if element is cut = enriched
           if cutlist(BOUNDARY(1,i).ele) == 0
              % element is not intersected by an interface

              % get nodenumbers of current boundary element
              nodenumbers = BOUNDARY(1,i).nodes;

              % get DOFs of the two nodes
              DOFs1 = id_eqns(BOUNDARY(i).nodes(1),1:2);
              DOFs2 = id_eqns(BOUNDARY(i).nodes(2),1:2);

              % get element load vector and ID array for
              % assembling into 'big_force'
              [force_values,force_id] = ...
                NBCintegrate(BOUNDARY(1,i),FORCE(1,j),DOFs1,DOFs2);

              % assemble into 'big_force'
              for k=1:4
                big_force(force_id(k)) = ...
                  big_force(force_id(k)) + force_values(k);
              end;
           else
%             if cutlist(BOUNDARY(1,i).ele) == 2
%               % element is intersected by one interface
% 
%               % get interface, that cuts the element
%               for k = 1:size(seg_cut_info,1)    % loop over interfaces 'k'
%                 for e = 1:size(seg_cut_info,2)  % loop over elements 'e'
%                   % if element 'e' is a boundary element
%                   if seg_cut_info(k,e).elemno == BOUNDARY(1,i).ele
%                     seg_cut_ind = [k e];
%                   end;
%                 end;
%               end;
% 
%               % get nodenumbers of current boundary element
%               nodenumbers = BOUNDARY(1,i).nodes;
% 
%               % get DOFs of the two nodes
%               DOFs1 = id_eqns(BOUNDARY(i).nodes(1),1:4);
%               DOFs2 = id_eqns(BOUNDARY(i).nodes(2),1:4);
% 
%               % get element load vector and ID array for
%               % assembling into 'big_force'
%               [force_values,force_id] = ...
%                 NBCintegrate_enriched(BOUNDARY(1,i),FORCE(1,j),DOFs1, ...
%                 DOFs2,seg_cut_info(seg_cut_ind(1),seg_cut_ind(2)), ...
%                 NODEINFO_ARR(BOUNDARY(1,i).nodes), ...
%                 id_dof(BOUNDARY(1,i).nodes,:));
% 
%               % assemble into 'big_force'
%               for k=1:2;%length(force_id)
%                 big_force(force_id(k)) = ...
%                   big_force(force_id(k)) + force_values(k);
%               end;
% 
%             end;
          end;
        end;
      end;
    end;
    otherwise
end;

% store maximum external load vector
big_force_max = big_force;

% clear some temporary variables
clear slot n j;
% ----------------------------------------------------------------------- %
%% INITIALZE (2)
% sum of all incremental solution vectors
fdisp_sum = zeros(size(big_force'));
% ----------------------------------------------------------------------- %
%% LOAD STEPPING LOOP (BEGIN)
% The deformed state will be computed via a incremental loading procedure.
% So, several interim states are computed. To be able to compare with the
% initial state, the inital x- and y-coordinates have to be stored in
% seperate variables. The nodal coordinates will be updated at the end of
% each load step iteration. The displacement increments will be saved to
% 'dis_increment'
x_orig = x; % x-coordinate
y_orig = y; % y-coordinate

% get a copy of the initial state of 'seg_cut_info'
seg_cut_info_orig = seg_cut_info;

% loop over all load / time steps (pseudo-time). The vector with time steps
% is given in the input file: 'IFtime'. Its first element is always a '0'.
for timestep = 1:(length(time)-1)
  disp(['load step ' num2str(timestep)]);
  % --------------------------------------------------------------------- %
  %% GET GLOBAL FORCE VECTOR FOR CURRENT LOAD STEP
  % get the current load increment for this loadstep, which is the basis to
  % assemble 'big_force' during the iterations. 'big_force' is built from
  % 'big_force_loadstep' (contribution of external loads = NBCs) and
  % 'big_force_traction' (contribution of internal tractions at the
  % interface due to plasticity or frictional sliding).
  big_force_loadstep = big_force_max * (time(timestep+1)-time(timestep));
  
  % This is also the global force vector for the first newton step.
  big_force = big_force_loadstep;
%   big_force(find(big_force))
  % --------------------------------------------------------------------- %
  %% SOLVE (QUASI-NEWTON-SCHEME BEGIN)
  % Due to the non-linearity of plasticity or friction, an iterative solver
  % is necessary. A Quasi-Newton-Scheme is applied.
  
  % The stiffness matrix is based on the elastic contribution 'bigk_el',
  % which is the same for all loadsteps, and some contributions due to the
  % interface, which depend on the sliding case 'IFsliding_switch'
  
  % 'fdisp' will be a vector with the solution for all global degrees of
  % freedom, with "base" and enriched degrees of freedom separated.

  % iterative solving via a Quasi-Newton-Scheme
  
  % initialize some variables
  solu = zeros(size(big_force));            % set a start vector
  iter = 0;                                 % iteration index for newton-scheme

  % store old solution to enable update
  old_solu = solu;
  
  while 1         % maximum number of iterations = maxiter
    % add '1' to the iteration index
    iter = iter + 1;
    % ------------------------------------------------------------------- %
    %% APPLY CONSTRAINTS AT INTERFACES (DEPENDING ON METHOD)

    % load elastic contribution of global stiffness matrix, which is 
    % identical in every loadstep
    bigk = bigk_el;
    
    % Now, add contributions due to the constraints at the interface. These
    % depend on the sliding case and the method of constraint enforcement.

    % Method is chosen by paramter 'IFmethod'
    switch IFmethod
      case 0  % Lagrange Multipliers
        % APPLY CONSTRAINTS AT GRAIN INTERFACES USING LAGRANGE MULTIPLIERS

%         disp('enforcing constraints at interfaces via Lagrange multipliers ...');

        % initialize
        ex_dofs = 0;
        lag_surf = [];
        
        count = 0;              % counter
        lag_slidestate = [];    % internal slidestate vector

        for i = 1:size(seg_cut_info,1)     % for every interface
          for e = 1:size(seg_cut_info,2) % for every cut element in that interface
            if seg_cut_info(i,e).elemno ~= -1
              % increase counter 'count'
              count = count + 1;
              
              % get global element ID
              parent_el = seg_cut_info(i,e).elemno;

              % Count the extra degree of freedom we're at
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
                gen_lagrange(node,x_orig,y_orig,parent_el,id_eqns,id_dof,pn_nodes, ...
                pos_g,neg_g,old_size,ex_dofs,seg_cut_info(i,e).xint, ...
                INTERFACE_MAP(i).endpoints);

              % stiffnes contribution 'ke_lag' has to be modified depending
              % on 'IFsliding_switch'
              switch sliding_switch
                case 0
                case 1              % frictionless sliding
                  % dot with normal vector
                  ke_lag = ke_lag * seg_cut_info(i,e).normal;
                  
                  % build full matrix
                  ke_lag_temp = [ke_lag zeros(12,1)];
                  
                  % assign full matrix
                  ke_lag = ke_lag_temp;
                  
                  % clear temporary variable
                  clear ke_lag_temp
                case 2              % perfect plasticity
                  % compute 'ke_pen' depending on current slide state (stick or slip)
                  if seg_cut_info(i,e).slidestate == 0      % stick
                    % equates to fully tied case --> no modification of
                    % elemetnary contribution to global stiffness matrix
                    
                    % set flag for internal slidestate-variable
                    lag_slidestate(count) = 0;
                  elseif seg_cut_info(i,e).slidestate == 1  % slip
                    % equates to frictionless sliding case
                                        
                    % dot with normal vector
                    ke_lag = ke_lag * seg_cut_info(i,e).normal;

                    % build full matrix
                    ke_lag_temp = [ke_lag zeros(12,1)];

                    % assign full matrix
                    ke_lag = ke_lag_temp;

                    % clear temporary variable
                    clear ke_lag_temp
                    
                    % set flag for internal slidestate-variable
                    lag_slidestate(count) = 1;
                  else
                    error('MATLAB:XFEM:UnvalidState', ...
                      'Current slide state is not valid!');
                  end;
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

        % Fix unused tangential equations ( = 'zero'-rows in 'bigk')
        switch IFsliding_switch 
          case 0  % fully tied case (no sliding at all)
            % In the fully tied case, there are no "unused tangential
            % eqautations" in 'bigk'.
          case 1  % frictionless sliding
            for i = 2:2:size(lag_surf,1)
              % get position index in 'bigk'
              global_index = i + max(max(id_eqns));
              
              % set '1' into 'bigk'
              bigk(global_index,global_index) = 1;
            end
          case 2  % perfect plasticity
            % reset counter
            count = 0;
            
            % loop over all constraint equations
            for i = 2:2:size(lag_surf,1)
              % increase counter
              count = count + 1;
              
              % set a '1' into 'bigk', if current constraint equation is a 
              % frictionless constraint (slip in current subsegment)
              if lag_slidestate(count) == 1
                % get position-index in 'bigk'
                global_index = i + max(max(id_eqns));
                
                % set '1' in 'bigk'
                bigk(global_index,global_index) = 1;
              end;
            end
          case 3              % frictional contact (Coulomb)
            warning('MATLAB:XFEM:main_xfem',...
              'There exists no code for frictional contact (Coulomb), yet.')
          otherwise
            warning('MATLAB:XFEM:main_xfem',...
              'Unvalid slidingID. Choose valid ID or add additional case to switch-case-structure')
        end;

        % clear temporary variable 'global_index'
        clear global_index i;
      case 1  % Penalty-Method
        % APPLY CONSTRAINTS AT GRAIN INTERFACES USING PENALTY METHOD
%         disp('enforcing constrainst at interfaces via penalty method ...');

        % get values from input file
        penalty = IFpenalty;

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
              %sliding cases are treated inside 'gen_penalty()'
              [ke_pen,id_pen] =...
                gen_penalty(node,x_orig,y_orig,parent_el,id_eqns,id_dof,...
                pn_nodes,pos_g,neg_g,seg_cut_info(i,e).xint,...
                INTERFACE_MAP(i).endpoints, ...
                seg_cut_info(i,e).normal,IFsliding_switch, ...
                seg_cut_info(i,e).slidestate);

              nlink = size(id_pen,2);
              %
              % assemble 'ke_pen' into 'bigk'
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
      case 2  % Nitsche's Method
        % APPLY CONSTRAINTS AT GRAIN INTERFACES USING NITSCHE'S METHOD
%         disp('enforcing constraints at interfaces via Nitsche´s method ...');

        % get values from input file
        penalty =IFnitsche;

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
              % treatment of different sliding cases inside 'gen_penalty()'
              [ke_pen,id_pen] =...
                gen_penalty(node,x_orig,y_orig,parent_el,id_eqns,id_dof,...
                pn_nodes,pos_g,neg_g,seg_cut_info(i,e).xint,...
                INTERFACE_MAP(i).endpoints, ...
                seg_cut_info(i,e).normal,IFsliding_switch, ...
                seg_cut_info(i,e).slidestate);                    

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
              % treatment of different sliding cases inside 'nit_stiff()'
              [ke_nit,id_nit] =...
                nit_stiff(node,x_orig,y_orig,parent_el,id_eqns,id_dof, ...
                pn_nodes,pos_g,neg_g,seg_cut_info(i,e).normal, ...
                seg_cut_info(i,e).xint,...
                INTERFACE_MAP(i).endpoints,IFsliding_switch,...
                seg_cut_info(i,e).slidestate);


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
    % ------------------------------------------------------------------- %
    %% APPLY DISPLACEMENT BOUNDARY CONDITIONS 
    % Due to an interative solver, prescribed displacements (unequal to 
    % zero) are only allowed in the first iteration step.
   
    for n=1:numnod  % loop over all nodes
      for j=1:ndof  % loop over all nodal DOFs
% The following commented block is the original one from Jessica Sanders,
% but it doesn't treat prescribed boundary conditions (unequal to zero) in
% combination with an iterative solver correctly.
%{
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
          elseif (m2 ~= 0) && (m3 == 0)             % If the node is enriched
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
%}
          % This is a modificated code, that skips the prescribed 
          % displacements for iteration steps iter > 1.
          if (dispbc(j,n) == 1)          
            m  = id_eqns(n,j);
            m2 = id_eqns(n,j+2);
            m3 = id_eqns(n,j+4);
            if (m2 == 0) && (m3 == 0)       % If the node is unenriched
              temp = size(bigk,2);
              bigk(m,:) = zeros(1,temp);
              if (iter == 1) || (ubar(j,n) == 0) && ...
                  ((IFsliding_switch == 0) || (IFsliding_switch == 1))
                big_force = big_force - (ubar(j,n) * bigk(:,m) * ...
                  (time(timestep + 1) - time(timestep)));
              elseif IFsliding_switch == 2
                big_force = big_force - (ubar(j,n) * bigk(:,m) * ...
                  (time(timestep + 1) - time(timestep)));
              end;
              
              bigk(:,m) = zeros(temp,1);
              bigk(m,m) = 1.0;
              big_force(m) = ubar(j,n) * (time(timestep + 1) - ...
                  time(timestep));
            elseif (m2 ~= 0) && (m3 == 0)             % If the node is enriched
              if id_dof(n,j+2) ~= nodegrainmap(n);    % But only 1 dof is active
                n;
                temp = size(bigk,2);
                bigk(m,:) = zeros(1,temp);
                if (iter == 1) || (ubar(j,n) == 0) && ...
                  ((IFsliding_switch == 0) || (IFsliding_switch == 1))
                  big_force = big_force - (ubar(j,n) * bigk(:,m) * ...
                    (time(timestep + 1) - time(timestep)));
                elseif IFsliding_switch == 2
                  big_force = big_force - (ubar(j,n) * bigk(:,m) * ...
                    (time(timestep + 1) - time(timestep)));
                end;
                
                bigk(:,m) = zeros(temp,1);
                bigk(m,m) = 1.0;
                big_force(m) = ubar(j,n) * (time(timestep+1) - ...
                    time(timestep));
              end;
            end;
          end;
      end;
    end;

    % clear some temporary variables
    clear m m2 m3 temp;

    % ------------------------------------------------------------------- %
    %% ADD CONSTRAINT EQUATIONS FOR DIRICHLET BOUNDARY CONDITIONS ON ENRICHED NODES
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
    % ------------------------------------------------------------------- %
    %% FIX NONPHYSICAL NODES (due to gmsh-meshes)
    % When using GMSH for mesh generation (unstructured meshes), there might be
    % some nodes, that are only needed for geometry purposes, not for meshing
    % (e.g. the center point of circles). These are referred to as "nonphysical
    % nodes". The corresponding equations have to be eliminated in the global
    % equations system.
    for i = nonphysnodevec  % 'nonphysnodevec' is empty, if there are no 
                            % non-physical nodes
      dofvec_temp = id_eqns(i,:); % get global DOF-numbers for nonphysical node
      for j=dofvec_temp
        bigk(j,j)=1;
      end;
    end;

    % clear some temporary variables
    clear dofvec_temp i j;
    % ------------------------------------------------------------------- %
    %% SOLVE (QUASI-NEWTON-SCHEME)
    % compute the increment vector 'delta'
    delta = -old_solu + bigk\big_force; 

    % update the solution-vector
    solu = solu + delta;

    % print some information about the current iteration step
    step_info = ['  Newton step: ' num2str(iter) '    Res-Norm: '...
        num2str(norm(delta))];
    disp(step_info);
    % ------------------------------------------------------------------- %
    %% CHECK YIELD CONDITIONS
    % A return mapping is not nevessary for the fully tied case and the
    % frictionless sliding case, but for all other cases.
    switch IFsliding_switch
      case 0  % fully tied case
        % The fully tied case is posed as a fully elastic problem with a
        % perfect bond. There are no yield conditions to check.
      case 1  % frictionless sliding
        % The frictionless sliding case is posed as a fully elastic problem
        % with a prefect bond in normal direction at the interface. There
        % are no yield conditions to check.
      case 2  % perfect plasticity
        % The plastic effects can only occur in tangential direction at the
        % interface. So, here a yield condition has to be checked to
        % decide, if a subsegment is in a 'stick' or a 'slip' state.
        %% COMPUTE CURRENT GLOBAL DISPLACEMENT VECTOR
        % compute global displacement vector as difference between current
        % and initial state
        for i=1:numnod                      % loop over all nodes 
          current_dis(2*i-1) = x(i) - x_orig(i);    % x-coordinate
          current_dis(2*i) = y(i) - y_orig(i);      % y-coordinate
        end;
        % --------------------------------------------------------------- %
        %% COMPUTE STRESSES ALONG THE INTERFACE
%         disp('compute stresses along the interface ...');

        % compute stresses 'stress_interface' and strains 'strain_interface' at
        % center of each intersected element along the interface. Stresses for
        % not-intersected elements are not computed, since they are not of
        % interest here.
        % Structure of 'stress_interface':
        %   Dimension: numelex6x3
        %   Index 1:    'stress_intercace' for element 'e'
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
        for i=1:size(seg_cut_info,1)            % loop over interface
          for e=1:size(seg_cut_info,2)          % loop over elements
            if seg_cut_info(i,e).elemno ~= -1   % only cut elements
              % get current element ID
              eleID = seg_cut_info(i,e).elemno;
              
              % compute stress in current element
              [straine,stresse] = post_process_better(node,x_orig, ...
                y_orig,eleID,current_dis,old_ndisp,id_dof,cutlist, ...
                maxngrains,INT_INTERFACE);
              
              % assign 'stresse' to 'stress_interface'
              stress_interface(eleID,1:6,:) = stresse;
              
              % assign 'straine' to 'strain_interface'
              strain_interface(eleID,1:6,:) = straine;
            end;
          end;
        end;

        % clear some temporary variables
        clear stresse straine i e;
        % --------------------------------------------------------------- %
        %% UPDATE SLIDING STATE AND TRACTION OF EACH INTERSECTED ELEMENT 
        % reset additional tractions at interface
        big_force_traction = zeros(length(big_force),1);
        
        % reset counter 'count'
        count = 0;
        
        % initialize 'slidestateconv' only in first newton-step
        if iter == 1
          % These flags indicate, if there are changes in the 'active set' 
          % in comparison with the previous iteration step 
          %   0 ... changes
          %   1 ... no changes
          for i=1:size(seg_cut_info,1)            % loop over interface
            for e=1:size(seg_cut_info,2)          % loop over elements
              if seg_cut_info(i,e).elemno ~= -1   % only cut elements
                % counter for intersected elements
                count = count + 1;

                % initialize to 'no changes'
                slidestateconv(count) = 1;  
                      % This will be checked later and changed, if necessary.
              end;
            end;
          end;
        end;
        
        % reset counter 'count'
        count = 0;
        
        % loop over all intersected elements and update the
        % 'slidestate'-flags
        for i=1:size(seg_cut_info,1)            % loop over interface
          for e=1:size(seg_cut_info,2)          % loop over elements
            if seg_cut_info(i,e).elemno ~= -1   % only cut elements
              % Now, get the maximum traction vector and the current
              % traction vector. Both are based on stress, since the yield
              % criterion is a stress criterion. Since the stresses on both
              % sides of the interface can be different (especially when
              % the materials are different), the greater traction vector
              % has to be used for further computation.
              
              % counter for intersected elements
              count = count + 1;
              
              % get current element ID
              eleID = seg_cut_info(i,e).elemno; 
              
              % get nodes of current element
              elenodes = node(:,eleID);

              % get global DOFs for 'elenodes'
              DOFs = id_eqns(elenodes,:);

              % get nodal contributions to 'big_force_traction' and the
              % vector with the maximum tangential traction (limited due to
              % plasticity)
              [force_values force_id tang_traction_max he] = ...
                getnodaltractionsplasticity(x_orig(elenodes), ...
                y_orig(elenodes),seg_cut_info(i,e),IFyieldstress, ...
                INTERFACE_MAP(i).endpoints,id_dof(elenodes,:),DOFs, ...
                NODEINFO_ARR(1,elenodes));
                   
              % get the two grain IDs of the current interface 'i'
              grains = seg_cut_info(i,e).grains;
                                 
              % compute current traction at interface depending on
              % 'Ifmethod'
              switch IFmethod
                case 0  % Lagrange multipliers 
                  % When a Lagrange multiplier method is used, the
                  % tractions at the interface are primary unknows and so
                  % they are part of the solution vector 'solu' or 'fdis'.
                  % To get the traction vector in a subsegment, the
                  % corresponding Lagrange multiplier has to be extracted
                  % form the solution vector.
                  
                  % Extract Lagrange multipliers from current
                  % solution-vector 'solu'
                  lagmult = fdisp_sum(old_size+1:old_size + 2 * multipliers) ...
                    + solu(old_size+1:old_size + 2 * multipliers)';

                  % find row in 'lagmult', that suits to current interface 'e'
                  % and to the element, that is cut by 'i'
                  index = find(lag_surf(:,2)==i & lag_surf(:,3)==eleID);
                  traction = lagmult(index);
                   
                  % reference frame of 'traction' depends on the slidestate
                  if seg_cut_info(i,e).slidestate == 0      % stick
                    % project traction into tangential direction  
                    tang_traction = (traction * ...
                      seg_cut_info(i,e).tangent) * seg_cut_info(i,e).tangent;
                  elseif seg_cut_info(i,e).slidestate == 1  % slip
                    tang_traction = traction';
                  else
                    error('MATLAB:XFEM:UnvalidState', ...
                      'Current slide state is not valid!');
                  end;
                                    
                  % clear some temporary variables
                  clear traction index lagmult;
                case 1  % penalty method
                  % Compute Lagrange multipliers via 
                  % 'lambda = alpha * \int _{\Gamma _e}[[u]] d \Gamma _e'

                  % get penalty parameter
                  penalty = IFpenalty;

                  % Establish which nodes are "postively" enriched, and
                  % which reside in the "negative" grain
                  pos_g = seg_cut_info(i,e).positive_grain;
                  neg_g = seg_cut_info(i,e).negative_grain;

                  [pn_nodes] =... 
                    get_positive_new(eleID,pos_g,neg_g);

                  % get traction vector at interface
                  traction = ...
                    get_lag_mults_for_penalty(node,x_orig,y_orig,eleID, ...
                    id_eqns,id_dof,pn_nodes,pos_g,neg_g, ...
                    seg_cut_info(i,e).xint, ...
                    INTERFACE_MAP(i).endpoints,penalty,fdisp_sum);

                  % project traction into tangential direction
                  tang_traction = (traction * ...
                    seg_cut_info(i,e).tangent) * seg_cut_info(i,e).tangent;
                  
%                   disp(num2str(norm(tang_traction)));

                  % clear some temporary variables
                  clear traction penalty pn_nodes pos_g neg_g;
                case 2  % Nitsche's method
                  % Compute Lagrange multipliers via 
                  % 'lambda = alpha * \int _{\Gamma _e}[[u]] d \Gamma _e ...
                  %                     - <sigma>*normal'

                  % get stabilization parameter
                  penalty = IFnitsche;
                  
                  grains = seg_cut_info(i,e).grains;
                  
                  % build stress tensors for current element
                  % in first grain
                  sigma_1 = [stress_interface(eleID,4,grains(1)) stress_interface(eleID,6,grains(1));
                    stress_interface(eleID,6,grains(1)) stress_interface(eleID,5,grains(1))];
                  % in second grain
                  sigma_2 = [stress_interface(eleID,4,grains(2)) stress_interface(eleID,6,grains(2));
                    stress_interface(eleID,6,grains(2)) stress_interface(eleID,5,grains(2))];
                  
                  % compute average stress
                  sigma_avg = (sigma_1 + sigma_2) / 2;
                  
                  % get nitsche contribution to traction vector at 
                  % interface
                  traction_nit = sigma_avg * seg_cut_info(i,e).normal;

                  % Establish which nodes are "postively" enriched, and
                  % which reside in the "negative" grain
                  pos_g = seg_cut_info(i,e).positive_grain;
                  neg_g = seg_cut_info(i,e).negative_grain;

                  [pn_nodes] =... 
                    get_positive_new(eleID,pos_g,neg_g);

                  % get penalty contribution to traction vector at 
                  % interface
                  traction_pen = ...
                    get_lag_mults_for_penalty(node,x_orig,y_orig,eleID, ...
                    id_eqns,id_dof,pn_nodes,pos_g,neg_g, ...
                    seg_cut_info(i,e).xint, ...
                    INTERFACE_MAP(i).endpoints,penalty,fdisp_sum);
                  
                  % compute traction vector at interface
                  traction = traction_pen - traction_nit';
                  
                  % project traction into tangential direction
                  tang_traction = (traction * ...
                    seg_cut_info(i,e).tangent) * seg_cut_info(i,e).tangent;
                
%                   disp(num2str(norm(tang_traction)));
                  
                  % clear some temporary variables
                    clear traction penalty pn_nodes pos_g neg_g ...
                      traction_pen traction_nit;
                otherwise
                  error('MATLAB:XFEM:UnvalidID',...
                    'Unvalid method ID. Choose valid ID or add additional case to switch-case-structure');
              end;
              
              % Determine sign of traction:
              % The tangential traction hast to be applied in the opposite 
              % direction of the tangential movement. Therefor, the last
              % displacement increment of the current element has to be
              % evaluated.
              
              % get displacement increment of last iteration step, averaged
              % at the center of the subsegment
              
              
              % set sign of 'force_values'
              if tang_traction' * tang_traction_max < 0
%                 disp('*********************************');
                force_values = force_values * (-1); 
              end;
              
              % check yield condition
              if norm(tang_traction) > norm(tang_traction_max)
                % "trial traction" > "yield traction" --> plastic
                % deformation
                %% check, if 'slidestate' changes
                if seg_cut_info(i,e).slidestate == 0;
                  slidestateconv(count) = 0;  % change in 'slidestate'
                else
                  slidestateconv(count) = 1;  % no change in 'slidestate'
                end;
                % ------------------------------------------------------- %
                %% set 'slidestate'
                seg_cut_info(i,e).slidestate = 1;
                % ------------------------------------------------------- %
                %% assemble new tangential traction into local force vector
                % 'big_force_traction'
                % assemble 'force_values' into 'big_force_traction'
                for j=1:length(force_values)
                  if force_id(j) ~= 0
                    big_force_traction(force_id(j)) = ...
                      big_force_traction(force_id(j)) + force_values(j);
                  end;
                end;
                % ------------------------------------------------------- %
              else
                % "trial traction" <= "yield traction" --> elastic
                % deformation
                %% check, if 'slidestate' changes
                if seg_cut_info(i,e).slidestate == 1;
                  slidestateconv(count) = 0;  % change in 'slidestate'
                else
                  slidestateconv(count) = 1;  % no change in 'slidestate' 
                end;
                % ------------------------------------------------------- %
                %% set 'slidestate'
                seg_cut_info(i,e).slidestate = 0;
                % ------------------------------------------------------- %
              end;
            end;
          end;
        end;
        
        % reset counter
        count = 0;
        
        % The assembling of the plastic tractions depends on the method of
        % constraint enforcement. 
        % Using Lagrange multipliers, there are extra variables for these
        % tractions, namely the Lagrange multipliers. For penalty and
        % Nitsche's method, where the displacements are the only  unknowns,
        % the tractions have to be assembled into the global force vector
        switch IFmethod
          case 0  % Lagrange multipliers
%             % loop over all interfaces 'i'
%             for i=1:size(seg_cut_info,1)
%               % loop over all cut elements
%               for e=1:size(seg_cut_info,2)
%                 % only over cut elements
%                 if seg_cut_info(i,e).elemno ~= -1
%                   % increase counter
%                   count = count + 1;
%                   
%                   % assemble the maximum tangential traction only into those
%                   % Lagrange multipliers in 'fdisp', which belong to an element
%                   % with 'slidestate = 1 [slip]'.
%                   if seg_cut_info(i,e).slidestate == 1
%                     % set 'flag' to '1' or '-1' to indicate the direction
%                     % of sliding
%                     if tang_traction' * tang_traction_max > 0
%                       flag = -1;
%                     else
%                       flag = 1;
%                     end;
%                         
%                     % assemble absolute value of maximum traction into
%                     % second (= tangential) Lagrange multiplier of this
%                     % element
%                     fdisp(old_size + 2 * count) = flag * ...
%                       norm(tang_traction_max);
%                   end;
%                 end;
%               end;
%             end;
          case 1  % penalty method
            % assemble interface tractions 'big_force_traction' into global 
            % force vector 'big_force'
            big_force = big_force_loadstep + big_force_traction * (time(timestep+1)-time(timestep));
          case 2  % Nitsche's method
            % assemble interface tractions 'big_force_traction' into global 
            % force vector 'big_force'
            big_force = big_force_loadstep + big_force_traction * (time(timestep+1)-time(timestep));
          otherwise
            error('MATLAB:XFEM:UnvalidID',...
              'Unvalid method ID. Choose valid ID or add additional case to switch-case-structure');
        end;
        
        % assemble interface tractions 'big_force_traction' into global 
        % force vector 'big_force'
        big_force = big_force_loadstep + big_force_traction * (time(timestep+1)-time(timestep));
        
        % clear some temporary variables
        clear tang_traction tang_traction_max stress_interface ...
          force_values force_id grains i e he eleID elenodes DOFs;
        % --------------------------------------------------------------- %
      case 3  % frictional contact (Coulomb)
        warning('MATLAB:XFEM:main_xfem',...
          'There exists no code for frictional contact (Coulomb), yet.')
      otherwise
        warning('MATLAB:XFEM:main_xfem',...
          'Unvalid slidingID. Choose valid ID or add additional case to switch-case-structure')
    end;
    % ------------------------------------------------------------------- %
    %% CHECK CONVERGENCE OF QUASI-NEWTON-SCHEME
    % show convergence state of 'active set'
%     disp(['slidestateconv = ' num2str(all(slidestateconv))]);
    
    % convergence check
    % The newton-scheme is assumed as converge, when the norm of the
    % displacement-increment is smaller than a given convergence tolerance
    % 'IFconvtol' and there are no changes in the sliding states in
    % comparison to the previous iteration step.
    if (norm(delta) < IFconvtol) && (all(slidestateconv))
      break;  % exit iteration loop, if convergence is achieved
    end;

    % check, if number of maximum iterations is reached
    if iter > IFmaxiter
      error('MATLAB:XFEM:main_xfem',...
          'Newton did not converge in %d iterations.', IFmaxiter);
    end;
    
    % store 'solu' to 'old_solu' for next iteration
    old_solu = solu;
    % ------------------------------------------------------------------- %
  end;
  % transpose 'solu'
  fdisp = solu';      % fdisp is a row-vector

  % sum 'fdisp' over all loadsteps
  fdisp_sum = fdisp_sum + fdisp;
  
  % clear some temporary variables
  clear solu step_info delta old_solu iter;
  % --------------------------------------------------------------------- %
  %% RE-ASSEMBLE GLOBAL DISPLACEMENT VECTOR
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
  
  % 'old_ndisp' stores displacement in base and enriched DOFs separated
  old_ndisp = old_ndisp + ndisp;    % update every load step
  
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
  % --------------------------------------------------------------------- %
  %% UPDATE AT END OF LOAD-STEP-LOOP
  % update nodal coordinates
  for i=1:numnod                % loop over all nodes
    x(i) = x(i) + dis(2*i-1);   % update x-coordinate
    y(i) = y(i) + dis(2*i);     % update y-coordinate
  end;
  
%   % update intersection points
%   for i=1:size(seg_cut_info,1)    % loop over all interfaces 'i'
%     for e=1:size(seg_cut_info,2)  % loop over all cut elements 'e'
%       if seg_cut_info(i,e).elemno ~= -1 % only, if there is no triple junction
%         % get global node IDs of nodes of current element
%         ele_nodes = node(:,seg_cut_info(i,e).elemno)
% 
%         % get initial x- and y-coordinates of nodes of the element
%         x_coords = x_orig(ele_nodes);
%         y_coords = y_orig(ele_nodes);
% 
%         % get global DOFs of these nodes
%         DOF_1 = id_eqns(node(1,seg_cut_info(i,e).elemno),1:2);  % first node
%         DOF_2 = id_eqns(node(2,seg_cut_info(i,e).elemno),1:2);  % second node
%         DOF_3 = id_eqns(node(3,seg_cut_info(i,e).elemno),1:2);  % third node
% 
%         % get total displacement for these nodes
%         dis_nodes = [dis(DOF_1); dis(DOF_2); dis(DOF_3)];
% 
%         % compute new intersection points
%         seg_cut_info(i,e).xint = ...
%           updateintersectionpoints(seg_cut_info_orig(i,e),x_coords, ...
%           y_coords,dis_nodes);
%       end;
%     end;
%   end;
  
  % store increment of displacement for each load step
  dis_increment(:,1,timestep) = dis;
  
  % updates at the end of each load step, that are specific with respect to
  % the sliding case
  switch IFsliding_switch
    case 0  % fully tied case
      % no updates necessary
    case 1  % frictionless sliding
      % no updates necessary
    case 2  % perfect plasticity
      % reset additional tractions at interface
      big_force_traction = zeros(length(big_force),1);
    case 3  % frictional contact (Coulomb)
      warning('MATLAB:XFEM:main_xfem',...
        'There exists no code for frictional contact (Coulomb), yet.')
    otherwise
      warning('MATLAB:XFEM:main_xfem',...
        'Unvalid slidingID. Choose valid ID or add additional case to switch-case-structure')
  end;
  % --------------------------------------------------------------------- %
  %% LOAD STEPPING LOOP (END)
end;
% ----------------------------------------------------------------------- %
%% COMPUTE GLOBAL DISPLACEMENT VECTOR
% compute global displacement vector as difference between deformed and
% initial state
for i=1:numnod                      % loop over all nodes 
  dis(2*i-1) = x(i) - x_orig(i);    % displacements of x-coordinates
  dis(2*i) = y(i) - y_orig(i);      % displacements of y-coordinates
end;
% ----------------------------------------------------------------------- %
%% LOAD INITIAL STATE
% load initial state to do some postprocessing.
x = x_orig; % x-coordinates of nodes
y = y_orig; % y-coordinates of nodes
% Now, the deformed state can be obtained by adding the displacement 'dis' 
% to the initial configuration 'x=x_orig', 'y=y_orig'.
% ----------------------------------------------------------------------- %
%% POST-PROCESS: STRESSES
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
        
  [straine,stresse] = post_process_better(node,x,y,e,dis,old_ndisp, ...
    id_dof,cutlist,maxngrains,INT_INTERFACE);
  stress(e,1:6,:) = stresse;
  strain(e,1:6,:) = straine;
end

% clear some temporary variables
clear stresse straine maxstress_vec minstress_vec i e f j;
% ----------------------------------------------------------------------- %
%% POST PROCESS: TRACTIONS AT INTERFACE
% computing the internal forces in the interface depends on the method of
% constraint enforcing
switch IFmethod 
  case 0  % Lagrange multipliers
    % extract vector with lagrange multipliers from 'fdisp'
    lagmult = fdisp_sum(old_size+1:old_size + 2*multipliers);

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
  case 1  % Penalty method
    % Compute Lagrange multipliers via 'lambda = alpha * [[u]]'

    % get penalty parameter
    penalty = IFpenalty;

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
             INTERFACE_MAP(i).endpoints,penalty,fdisp_sum);

          lagmult = [lagmult seg_cut_info(i,e).lagmult];
        end;
      end;
    end;

    % clear some temporary variables
    clear pn_nodes neg_g pos_g i e parent_el;
  case 2  % Nitsche's method
    % Compute Lagrange multipliers via 
    %                       'lambda = alpha * [[u]] - <sigma>*normal'

    % get stabilization parameter
    penalty = IFnitsche;

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
             
          % grains, associated to interface 'i'
          grain1 = seg_cut_info(i,e).grains(1);
          grain2 = seg_cut_info(i,e).grains(2);

          % stress tensors in cut element 'e' on both sides of
          % interface 'i'
          sigma1 = [stress(parent_el,4,grain1) stress(parent_el,6,grain1);
            stress(parent_el,6,grain1) stress(parent_el,5,grain1)];
          sigma2 = [stress(parent_el,4,grain2) stress(parent_el,6,grain2);
            stress(parent_el,6,grain2) stress(parent_el,5,grain2)];

          % compute Nitsche Lagrange multiplier
          sigma_avg = 0.5 * (sigma1 + sigma2);
                    
          % get normal to the interface
          normal = seg_cut_info(i,e).normal;
          
          % compute penalty part of Lagrange multiplier as in
          % penalty case
          lagmult_pen = ...
            get_lag_mults_for_penalty(node,x,y,parent_el, ...
            id_eqns,id_dof,pn_nodes,pos_g,neg_g, ...
            seg_cut_info(i,e).xint, ...
            INTERFACE_MAP(i).endpoints,penalty,fdisp_sum);
              
          % compute Nitsche-part
          lagmult_nit = sigma_avg * normal;
                        
          % subtract the nitsche form the penalty contribution 
          seg_cut_info(i,e).lagmult = lagmult_pen - lagmult_nit';
                    
          % store into global vector 'lagmult' to find maximum
          % and minimum value afterwards
          lagmult = [lagmult seg_cut_info(i,e).lagmult];
        end;
      end;
    end;

    % clear some temporary variables
    clear pn_nodes neg_g pos_g i e parent_el lagmult_nit lagmult_pen ...
      normal sigma_avg sigma1 sigma2 grain1 grain2;
  otherwise
    error('MATLAB:XFEM:UnvalidID',...
      'Unvalid method ID. Choose a valid ID or add an additional case to switch-case-structure.');
end;
% ----------------------------------------------------------------------- %
%% FINISH SOLVING PROCESS
%  disp('saving to results file ...');

%  save my_results_files.mat x y node dis fdisp numele numnod cutlist...
%      INT_INTERFACE SUBELEM_INFO SUBELEMENT_GRAIN_MAP id_eqns id_dof...
%      X Y CONN stress

disp('Solving and postprocessing done.');

