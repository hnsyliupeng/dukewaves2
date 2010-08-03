%% Two dimensional polycrystal XFEM code
% G-FEM Based on the paper by Simone, Duarte and Van der Giessen (2006)
% All other ideas by Dolbow
% Jessica Sanders, Summer 2006 (SNL) - Summer 2007
% Furhter developements by Matthias Mayr (Summer 2010)
%% LOAD INPUT PARAMETERS FROM 'xfeminputdata_xfem.mat'
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
if exist('IFintegral','var') == 0, IFintegral = 2;end;      % variante of penalty term

if exist('IFyieldstress','var') == 0, IFyieldstress = 0;end;% yield stress

if exist('IFpenalty_normal','var') == 0, IFpenalty_normal = IFpenalty;end;  % Penalty-Parameter for normal direction
if exist('IFpenalty_tangential','var') == 0, IFpenalty_tangential = IFpenalty_normal;end; % Penalty-Parameter for tangential direction

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
%   dis_increment = zeros(2*numnod,1,(length(time)-1));
else
  time = [0 1];
%   dis_increment = zeros(2*numnod,1,1);
end;

% load pseudo-time-vectors, if additional sets of DBCs are defined
if exist('IFtime2','var')
  time2 = IFtime2;
end;
if exist('IFtime3','var')
  time3 = IFtime3;
end;
% ----------------------------------------------------------------------- %
%% INITIALIZE (1)
% counter 
count = 0;

% vector for 'Lagrange multipliers' (= tractions at interface)
lagmult = [];

% initialize a standard displacement vector
dis = zeros(1,2*numnod);

% initialize some variables which are needen in every subsegment of the
% interface
%   initialize sliding state flags 'slidestate' in 'seg_cut_info'
%     flag    description
%     0       stick
%     1       slip
%   initialize a variable for the flux via domain integral, too.
%   current tangential plastic gap 'tgappl' at the 2 gauss points
%   tangential plastic gap of previous laod step 'tgapplconv' at the 2
%     gauss points
%   final tangential gap of the previous load step 'tgapconv'
%   current tangential traction 'ttrac' at the 2 gauss points
%   final tangential traction of previous load step 'ttracconv' at the two
%     gauss points
%   evalutated flow rules of trial state 'f_trial' at both gauss points
for i=1:size(seg_cut_info,1)
  for e=1:size(seg_cut_info,2)
    seg_cut_info(i,e).slidestate = 0;     % initialize as 'stick'
%     seg_cut_info(i,e).domint = [];        % vector for traction at interface
    seg_cut_info(i,e).tgappl = [0 0];     % current plastic gap
    seg_cut_info(i,e).tgapplconv = [0 0]; % converged plastic gap of previous load step
%     seg_cut_info(i,e).tgapconv = [0 0];   % converged total gap of previous load step
    seg_cut_info(i,e).ttrac = [0 0];  % current tangential traction    
    seg_cut_info(i,e).ttracconv = [0 0];  % converged tangential traction of previous load step
    seg_cut_info(i,e).f_trial = [0 0];    % evalutated flow rules of trial state at both gauss points
  end;
end;

% initialize some variables related to plasticity which are needed to make 
% the code consistent for pure elastic problems
% initialize tangential plastic gap
tgappl = [0 0];

% % initialize a state-variable for normal contact
% % 0 ... contact open (tension, gap > 0)
% % 1 ... contact closed (pressure), gap < 0)
% for i=1:size(seg_cut_info,1)
%   for e=1:size(seg_cut_info,2)
%     seg_cut_info(i,e).contactstate = 1; % initialize as 'closed'
%   end;
% end;

% force-vectors for domain integral method
domint_force = zeros(18,1);

% stiffness-matrix for domain integral method
domint_stiff = zeros(18,18);

% displacement-vector for domain integral method
domint_dis = zeros(18,1);
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
  %% initialize
  % 'old_size' is the number of base-DOFs plus number of enriched DOFs 
  old_size = numeqns;

  % initialize the total number of extra added equations
  multipliers = 0;
  % --------------------------------------------------------------------- %
  %% add extra equations (only for Lagrange multipliers)
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
  % --------------------------------------------------------------------- %
  %% add extra equations to enforce DBCs on enriched nodes via Lagrange
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
  % --------------------------------------------------------------------- %
  %% allocate size of the stiffness matrix
  %bigk = zeros(numeqns);
  bigk = spalloc(numeqns,numeqns,numele*36);
  ndof = 2; % number of standard degrees of freedon per node
  
  % clear some temporary variables
  clear i e tempnode;
  % --------------------------------------------------------------------- %
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
clear rbk cbk re ce nlink id ke i e j bigk;
% ----------------------------------------------------------------------- %
%% PREPARE ASSEMBLY OF EXTERNAL FORCE VECTOR
disp('assembling external force vector ...')
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
% big_force_traction = zeros(numeqns,1);
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
            if cutlist(BOUNDARY(1,i).ele) == 2
              % element is intersected by one interface

              % get interface, that cuts the element
              for k = 1:size(seg_cut_info,1)    % loop over interfaces 'k'
                for e = 1:size(seg_cut_info,2)  % loop over elements 'e'
                  % if element 'e' is a boundary element
                  if seg_cut_info(k,e).elemno == BOUNDARY(1,i).ele
                    seg_cut_ind = [k e];
                  end;
                end;
              end;

              % get nodenumbers of current boundary element
              nodenumbers = BOUNDARY(1,i).nodes;

              % get DOFs of the two nodes
              DOFs1 = id_eqns(BOUNDARY(i).nodes(1),1:4);
              DOFs2 = id_eqns(BOUNDARY(i).nodes(2),1:4);

              % get element load vector and ID array for
              % assembling into 'big_force'
              [force_values,force_id] = ...
                NBCintegrate_enriched(BOUNDARY(1,i),FORCE(1,j),DOFs1, ...
                DOFs2,seg_cut_info(seg_cut_ind(1),seg_cut_ind(2)), ...
                NODEINFO_ARR(BOUNDARY(1,i).nodes), ...
                id_dof(BOUNDARY(1,i).nodes,:));

              % assemble into 'big_force'
              for k=1:2;%length(force_id)
                big_force(force_id(k)) = ...
                  big_force(force_id(k)) + force_values(k);
              end;

            end;
          end;
        end;
      end;
    end;
  otherwise
end;

% store maximum external load vector
big_force_max = big_force;

% clear some temporary variables
clear slot n j big_force DOFs1 DOFs2 nodenumbers force_values force_id k e i;
% ----------------------------------------------------------------------- %
%% INITIALZE (2)
% total displacement vector of all DOFs in the inner of newton loop
totaldis = zeros(size(big_force_max));

% needed to compute stresses
orig_ndisp = zeros(numnod,6);
% ----------------------------------------------------------------------- %
%% LOAD STEPPING LOOP (BEGIN)
% The deformed state will be computed via a incremental loading procedure.
% So, several interim states are computed. 

% loop over all load / time steps (pseudo-time). The vector with time steps
% is given in the input file: 'IFtime'. Its first element is always a '0'.
for timestep = 1:(length(time)-1)
  % print number of current load step
  loadsteptext = sprintf('\nlaod step %d',timestep);
  disp(loadsteptext);
  % --------------------------------------------------------------------- %
  %% GET GLOBAL FORCE VECTOR FOR CURRENT LOAD STEP
  % get the current load increment for this loadstep, which is the basis to
  % assemble 'big_force' during the iterations. 'big_force' is built from
  % 'big_force_loadstep' (contribution of external loads = NBCs) and
  % 'big_force_traction' (contribution of internal tractions at the
  % interface due to plasticity or frictional sliding).
  big_force_loadstep = big_force_max * time(timestep+1);
  % --------------------------------------------------------------------- %
  %% SOLVE (QUASI-NEWTON-SCHEME BEGIN)
  % Due to the non-linearity of plasticity or friction, an iterative solver
  % is necessary. A Quasi-Newton-Scheme is applied.
  
  % The stiffness matrix is based on the elastic contribution 'bigk_el',
  % which is the same for all loadsteps, and some contributions due to the
  % interface, which depend on the sliding case 'IFsliding_switch'
  
  % 'fdisp = totaldis' ' will be a vector with the solution for all global degrees of
  % freedom, with "base" and "enriched" degrees of freedom separated.

  % iterative solving via a Quasi-Newton-Scheme
  
  % initialize some variables
  iter = 1;                         % iteration index for newton-scheme
  slidestateconv = 'false';         % shows, if the 'slidestate'-flags 
                                    % are converged
%   contactstateconv = 'false';       % shows, if the 'contactstate'-flags 
                                    % are converged
  deltanewton = zeros(size(big_force_max));   % displacement increment for Newton loop
  deltaload = zeros(size(big_force_max));     % displacement increment for load step loop
  
  % initialize 'slidestateconv' only in first Newton-step
  % These flags indicate, if there are changes in the 'active set' 
  % in comparison with the previous iteration step 
  %   0 ... changes
  %   1 ... no changes
  for i=1:size(seg_cut_info,1)            % loop over interface
    for e=1:size(seg_cut_info,2)          % loop over elements
      if seg_cut_info(i,e).elemno ~= -1   % only cut elements
        % initialize
        seg_cut_info(i,e).slidestateconv = 0;  
%         seg_cut_info(i,e).contactstateconv = 0; 
        
%         if timestep == 1
          % initialize all elements to 'stick'
          seg_cut_info(i,e).slidestate = 0;
%         end;

          % initialize all elements to a closed contact
%           seg_cut_info(i,e).contactstate = 1;
      end;
    end;
  end;
    
  % loop for Quasi-Newton-Scheme
  while 1%(res_norm > IFconvtol && delta_norm > IFconvtol && ...
  %     ener_norm > IFconvtol) || strcmp(slidestateconv,'false') || ...
  %     strcmp(contactstateconv,'false')
  
%   % reset slidestate-flag
%   for i=1:size(seg_cut_info,1)            % loop over interface
%     for e=1:size(seg_cut_info,2)          % loop over elements
%       seg_cut_info(i,e).slidestate = 0;
%     end;
%   end;
  for i=1:size(seg_cut_info,1)            % loop over interface
    for e=1:size(seg_cut_info,2)          % loop over elements
      seg_cut_info(i,e).f_trial = [0 0];
    end;
  end;
%   % --------------------------------------------------------------------- %
    %% PREPARE IMPOSING DIRICHLET BOUNDARY CONDITIONS
    % The global displacement-vector 'totaldis' is copied and modified such
    % that its copy contains all Dirichlet boundary conditions for the
    % current load step. This modified displacement vector will be used to
    % assemble the residual.
 
    % initialize a vector to store all prescribed displacements
    DBCmatrix = zeros(2,length(totaldis));  % first line shows, if there is a DBC
                                            % seocnd line stores the value
    
    % loop over all nodes
    for n=1:numnod
      % loop x- and y-DOF of each node
      for j=1:2
        m  = id_eqns(n,j);    % global DOF ID for base DOF
%         m2 = id_eqns(n,j+2);  % global DOF ID for first enriched DOF
%         if m2 ~= 0  
%           m,m2,n,j% check, if node is enriched
%           error('MATLAB:XFEM:DBCs', ...
%             'No Dirichlet boundary conditions allowed on enriched nodes.');
%         end;
        
        % build vector with prescribed displacements
        if (dispbc(j,n) == 1)
          DBCmatrix(1,m) = 1;
          DBCmatrix(2,m) = ubar(j,n) * time(timestep + 1);
        end
        
        % add a second set of DBCs, if defined
        if exist('time2','var') == 1
          if (dispbc2(j,n) == 1)
            DBCmatrix(1,m) = 1;
            DBCmatrix(2,m) = DBCmatrix(2,m) + ubar2(j,n) * time2(timestep + 1);
          end;
        end;
        
        % add a third set of DBCs, if defined
        if exist('time3','var') == 1
          if (dispbc3(j,n) == 1)
            DBCmatrix(1,m) = 1;
            DBCmatrix(2,m) = DBCmatrix(2,m) + ubar3(j,n) * time3(timestep + 1);
          end;
        end;
      end;
    end;
    
    % build the modified displacement vector 'totaldisDBC'
    totaldisDBC = totaldis;               % copy global displacement vector
    for m=1:length(totaldis)              % loop over all DOFs
      if DBCmatrix(1,m) == 1              % check, if DOF has a DBC
        totaldisDBC(m) = DBCmatrix(2,m);  % assign its value
      end;
    end;
    
    % clear some temporary variables
    clear n j m;
    % ------------------------------------------------------------------- %
    %% BUILD RESIDUAL
%     % The global stiffnes matrix is already built and the Dirichlet
%     % boundary conditions are imposed. A global force vector 'big_force' is
%     % assembled. So, now the residual 'residual' can be build, whereby the
%     % internal force is computed based on stiffness matrix 'bigk' and the
%     % current global total displacement vector 'totaldis'.
%     residual = big_force - bigk * totaldis;

    % initialize residual
    residual = zeros(length(big_force_max),1);

%{
    % Assemble the residual contribution of internal forces on an element 
    % based level.
    for e = 1:numele    % loop over all elements
      % get some element data
      elenodes = node(:,e);     % global node IDs of this element's nodes
      xcoords = x(elenodes);    % x-coordinates of 'elenodes'
      ycoords = y(elenodes);    % y-coordinates of 'elenodes'
      
      % get id-array to prepare assembly
      id = [id_eqns(elenodes(1),1:2) id_eqns(elenodes(2),1:2) id_eqns(elenodes(3),1:2)...
          id_eqns(elenodes(1),3:6) id_eqns(elenodes(2),3:6) id_eqns(elenodes(3),3:6)];
      eliminate = find(id == 0);
      for i = size(eliminate,2):-1:1  % eliminate dofs with index '0'
        id(eliminate(i)) = [];
      end
      
      % Compute contributions to an elements residual
%       ele_residual_int = get_ele_residual_int(elenodes,xcoords,ycoords, ...
%         dis,orig_ndisp,cutlist(e),maxngrains,e,id_dof(elenodes,:));          % residual of internal forces due to the inner of the grains
%       ele_residual_neumann = get_ele_residual_neumann(xcoords,ycoords);  % residual of tractions on Neumann boundary
%       ele_residual_body = residual_body();        % residual of body forces

      if cutlist(e) == 0
        nodal_dis = [ totaldis(id_eqns(elenodes(1),1)); ...
                      totaldis(id_eqns(elenodes(1),2)); ...
                      totaldis(id_eqns(elenodes(2),1)); ...
                      totaldis(id_eqns(elenodes(2),2)); ...
                      totaldis(id_eqns(elenodes(3),1)); ...
                      totaldis(id_eqns(elenodes(3),2))];
        ele_residual = elemstiff_class * nodal_dis;
      else
        nodal_dis = [ totaldis(id_eqns(elenodes(1),1)); ...
                      totaldis(id_eqns(elenodes(1),2)); ...
                      totaldis(id_eqns(elenodes(2),1)); ...
                      totaldis(id_eqns(elenodes(2),2)); ...
                      totaldis(id_eqns(elenodes(3),1)); ...
                      totaldis(id_eqns(elenodes(3),2))];
      end;
      % compute element's residual
      ele_residual = ele_residual_int;% - ele_residual_neumann;
           
      % assemble into global residual
      for index = 1:length(id)
        residual(id(index)) = residual(id(index)) + ele_residual(index);
      end;
    end; 
%}
    
    % contribution of elastic bulk field of the inner of the grains which
    % will not be affected from plasticity at the interface
    residual = bigk_el * totaldisDBC;
    
    % residual contributions due to constraints
    for i=1:size(seg_cut_info,1)    % loop over all interfaces 'i'
      for e=1:size(seg_cut_info,2)  % loop over all cut elements 'e'
        if seg_cut_info(i,e).elemno ~= -1  % only, if 'e' is cut by 'i'
          % get some element data
          eleID = seg_cut_info(i,e).elemno; % global element ID
          elenodes = node(:,eleID); % global node IDs of this element's nodes
          xcoords = x(elenodes);    % x-coordinates of 'elenodes'
          ycoords = y(elenodes);    % y-coordinates of 'elenodes'
          
          switch IFmethod
            case 0  % Lagrange multipliers
            case 1  % Penalty method
              % get penalty parameters (depending on sliding case)
              switch IFsliding_switch
                case 0  % fully tied case
                  penalty_normal = IFpenalty_normal;
                  penalty_tangent = IFpenalty_tangential;
                case 1  % frictionless sliding
                  penalty_normal = IFpenalty_normal;
                  penalty_tangent = 0;
                case 2  % perfect plasticity
                  % provide both penalty parameters: the normal one will be
                  % used indepentent of the slidestate, the tangential one
                  % is needed for the return mapping algorithm.
                  penalty_normal = IFpenalty_normal;
                  penalty_tangent = IFpenalty_tangential;
                otherwise
                  error('MATLAB:XFEM','Unvalid Sliding-ID');
              end;
              
              % get residual contributions vectors, but choose betweeb the
              % one- and two-integral formulation
              switch IFintegral
                case 1  % one integral
                  [res_penalty_normal res_penalty_tangent id tgappl ttrac f_trial] = ...
                    get_ele_residual_penalty_alternative(xcoords,ycoords, ...
                    seg_cut_info(i,e),INTERFACE_MAP(i).endpoints, ...
                    id_dof(elenodes,:),id_eqns(elenodes,:),totaldis, ...
                    penalty_normal,penalty_tangent,IFyieldstress, ...
                    deltaload,IFsliding_switch);
                case 2  % two integrals
                  [res_penalty_normal res_penalty_tangent id] = ...
                    get_ele_residual_penalty(xcoords,ycoords, ...
                    seg_cut_info(i,e),INTERFACE_MAP(i).endpoints, ...
                    id_dof(elenodes,:),id_eqns(elenodes,:),totaldis);
                otherwise
                  error('MATLAB:XFEM:UnvalidID', ...
                    'Unvalid number of integrals');
              end;

              % store plastic contribution to tangential gap
              seg_cut_info(i,e).tgappl = tgappl;
              
              % store scalar values of current tangential traction at each
              % gauss point
              seg_cut_info(i,e).ttrac = ttrac;
              
              % store evaluated flow rules of trial state at both gauss
              % points
              seg_cut_info(i,e).f_trial = f_trial;
              
              % compute ele_residual_constraint
              ele_residual_constraint = res_penalty_normal ...  % normal traction
                                      + res_penalty_tangent;    % tangential traction
                                                                 
              % assemble into global residual
              for a = 1:length(id)
                if id(a) ~= 0
                  residual(id(a)) = residual(id(a)) + ele_residual_constraint(a);
                end;
              end;
            case 2  % Nitsche's mehtod
            otherwise
              error('MATLAB:XFEM','Unvalid Method-ID');
          end;
        end;
      end;
    end;

    residual = residual - big_force_loadstep;
    
    % clear some temporary variables
    clear eleID elenodes xcoords ycoords res_constant_tangent i e a id ...
      ele_residual_constraint constant_tangent res_penalty_normal ...
      penalty_normal res_penalty_tangent penalty_tangent ...
      normaltraction tangentialtraction;
    % ------------------------------------------------------------------- %
    %% BUILD TANGENT MATRIX
    % Here, the tangent matrix which is computed via
    % Dresidual/Ddisplpacements is assembled. Since the material in the
    % grains is assumed as linear elastic, its linearization is independent
    % from the displacements due to the small perturbation setting. So, it
    % is computed before the Newton loop and only copied here. It is used
    % as the starting part for 'tangentmatrix'.
    tangentmatrix = bigk_el;  
    
    % Now, additional contributions have to be assembled. First, assemble
    % the contributions to enforce the constraints. Of course, this depends
    % on the method of constraint enforcement:
    switch IFmethod
      case 0  % Lagrange multipliers
        error('MATLAB:XFEM:UnvalidID', ...
          'No code for Lagrange multipliers, yet');
      case 1  % Penalty method
        % loop over all interfaces and cut elements and assemble their 
        % penalty contributions into 'tangentmatrix'
        for i=1:size(seg_cut_info,1)
          for e=1:size(seg_cut_info,2)
            if seg_cut_info(i,e).elemno ~= -1
              % get some element data
              eleID = seg_cut_info(i,e).elemno;   % global element ID
              elenodes = node(:,eleID);           % global node IDs
              xcoords = x(elenodes);              % x-coordinates
              ycoords = y(elenodes);              % y-coordinates
              
              % distinguish between sliding cases
              switch IFsliding_switch
                case 0  % fully tied case
                  penalty_normal = IFpenalty_normal;
                  penalty_tangent = IFpenalty_tangential;
                case 1  % frictionless sliding
                  penalty_normal = IFpenalty_normal;
                  penalty_tangent = 0;
                case 2  % perfect plasticity
                  penalty_normal = IFpenalty_normal;
                  penalty_tangent = IFpenalty_tangential;
                otherwise
                  error('MATLAB:XFEM:UnvalidID', ...
                    'Unvalid sliding ID');
              end;
                            
              % get penalty matrices, but choose between one- and
              % two-integral formulation
              switch IFintegral
                case 1  % one integral
                  [pen_normal pen_tangent id_pen] = ... 
                    lin_penalty_alternative(xcoords,ycoords, ...
                    id_eqns(elenodes,:),id_dof(elenodes,:), ...
                    seg_cut_info(i,e),INTERFACE_MAP(i).endpoints, ...
                    penalty_normal,penalty_tangent);
                case 2  % two integrals
                  [pen_normal pen_tangent id_pen] = ... 
                    lin_penalty(xcoords,ycoords,id_eqns(elenodes,:), ...
                    id_dof(elenodes,:),seg_cut_info(i,e), ...
                    INTERFACE_MAP(i).endpoints);
                otherwise
                  error('MATLAB:XFEM:UnvalidID', ...
                    'Unvalid number of integrals');
              end;
              
              % build an elemental matrix
              penalty_matrix  = pen_normal + pen_tangent;
                            
              % assemble 'penalty_matrix' into 'tangentmatrix'
              nlink = size(id_pen,2);
              for m=1:nlink
                for n=1:nlink
                  rbk = id_pen(m);
                  cbk = id_pen(n);
                  re = m;
                  ce = n;
                  if ((rbk ~= 0) && (cbk ~= 0))                        
                    tangentmatrix(rbk,cbk) = tangentmatrix(rbk,cbk) + ...
                        penalty_matrix(re,ce);
                  end;
                end;
              end;
            end;
          end;
        end;
      case 2  % Nitsche's method
        error('MATLAB:XFEM:UnvalidID', ...
          'No code for Nitsche´s method, yet');
      otherwise
        error('MATLBA:XFEM:UnvalidID','Unvalid method-ID');
    end;
    % ------------------------------------------------------------------- %
    %% IMPOSE DISPLACEMENT BOUNDARY CONDITIONS 
    % loop over all DOFs
    numdof = length(totaldis);
    for m=1:numdof
      % check, if current DOF has a DBC
      if DBCmatrix(1,m) == 1
        tangentmatrix(:,m) = zeros(numdof,1);
        tangentmatrix(m,:) = zeros(1,numdof);
        tangentmatrix(m,m) = 1.0;
        residual(m) = totaldis(m) - totaldisDBC(m);
      end;
    end;

    % clear some temporary variables
    clear m numdof DBCmatrix totaldisDBC;
    % ------------------------------------------------------------------- %
    %% ADD CONSTRAINT EQUATIONS FOR DIRICHLET BOUNDARY CONDITIONS ON ENRICHED NODES
    %{
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
    %}
    % ------------------------------------------------------------------- %
    %% FIX NONPHYSICAL NODES (due to gmsh-meshes)
    %{
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
    %}
    % ------------------------------------------------------------------- %
    %% SOLVE (NEWTON-RAPHSON-SCHEME)
    % compute the increment vector 'deltanewton'
    deltanewton = -tangentmatrix\residual; % no negative sign here, since 'residual' 
                                 % is considered as 'F_ext - F_int' 
                                 % (according to Laursen's Book)
    % ------------------------------------------------------------------- %
%     %% UPDATE RESIDUAL
% %     % The global stiffnes matrix is already built and the Dirichlet
% %     % boundary conditions are imposed. A global force vector 'big_force' is
% %     % assembled. So, now the residual 'residual' can be build, whereby the
% %     % internal force is computed based on stiffness matrix 'bigk' and the
% %     % current global total displacement vector 'totaldis'.
% %     residual = big_force - bigk * totaldis;
%     
%     % copy old residual
%     residual_old = residual;
% 
%     % initialize residual
%     residual = zeros(length(big_force_max),1);
% 
% %{
%     % Assemble the residual contribution of internal forces on an element 
%     % based level.
%     for e = 1:numele    % loop over all elements
%       % get some element data
%       elenodes = node(:,e);     % global node IDs of this element's nodes
%       xcoords = x(elenodes);    % x-coordinates of 'elenodes'
%       ycoords = y(elenodes);    % y-coordinates of 'elenodes'
%       
%       % get id-array to prepare assembly
%       id = [id_eqns(elenodes(1),1:2) id_eqns(elenodes(2),1:2) id_eqns(elenodes(3),1:2)...
%           id_eqns(elenodes(1),3:6) id_eqns(elenodes(2),3:6) id_eqns(elenodes(3),3:6)];
%       eliminate = find(id == 0);
%       for i = size(eliminate,2):-1:1  % eliminate dofs with index '0'
%         id(eliminate(i)) = [];
%       end
%       
%       % Compute contributions to an elements residual
% %       ele_residual_int = get_ele_residual_int(elenodes,xcoords,ycoords, ...
% %         dis,orig_ndisp,cutlist(e),maxngrains,e,id_dof(elenodes,:));          % residual of internal forces due to the inner of the grains
% %       ele_residual_neumann = get_ele_residual_neumann(xcoords,ycoords);  % residual of tractions on Neumann boundary
% %       ele_residual_body = residual_body();        % residual of body forces
% 
%       if cutlist(e) == 0
%         nodal_dis = [ totaldis(id_eqns(elenodes(1),1)); ...
%                       totaldis(id_eqns(elenodes(1),2)); ...
%                       totaldis(id_eqns(elenodes(2),1)); ...
%                       totaldis(id_eqns(elenodes(2),2)); ...
%                       totaldis(id_eqns(elenodes(3),1)); ...
%                       totaldis(id_eqns(elenodes(3),2))];
%         ele_residual = elemstiff_class * nodal_dis;
%       else
%         nodal_dis = [ totaldis(id_eqns(elenodes(1),1)); ...
%                       totaldis(id_eqns(elenodes(1),2)); ...
%                       totaldis(id_eqns(elenodes(2),1)); ...
%                       totaldis(id_eqns(elenodes(2),2)); ...
%                       totaldis(id_eqns(elenodes(3),1)); ...
%                       totaldis(id_eqns(elenodes(3),2))];
%       end;
%       % compute element's residual
%       ele_residual = ele_residual_int;% - ele_residual_neumann;
%            
%       % assemble into global residual
%       for index = 1:length(id)
%         residual(id(index)) = residual(id(index)) + ele_residual(index);
%       end;
%     end; 
% %}
%     
%     % contribution of elastic bulk field of the inner of the grains which
%     % will not be affected from plasticity at the interface
%     residual = bigk_el * totaldisDBC;
%     
%     % residual contributions due to constraints
%     for i=1:size(seg_cut_info,1)    % loop over all interfaces 'i'
%       for e=1:size(seg_cut_info,2)  % loop over all cut elements 'e'
%         if seg_cut_info(i,e).elemno ~= -1  % only, if 'e' is cut by 'i'
%           % get some element data
%           eleID = seg_cut_info(i,e).elemno; % global element ID
%           elenodes = node(:,eleID); % global node IDs of this element's nodes
%           xcoords = x(elenodes);    % x-coordinates of 'elenodes'
%           ycoords = y(elenodes);    % y-coordinates of 'elenodes'
%           
%           switch IFmethod
%             case 0  % Lagrange multipliers
%             case 1  % Penalty method
%               % get penalty parameters (depending on sliding case)
%               switch IFsliding_switch
%                 case 0  % fully tied case
%                   penalty_normal = IFpenalty_normal;
%                   penalty_tangent = IFpenalty_tangential;
%                 case 1  % frictionless sliding
%                   penalty_normal = IFpenalty_normal;
%                   penalty_tangent = 0;
%                 case 2  % perfect plasticity
%                   % provide both penalty parameters: the normal one will be
%                   % used indepentent of the slidestate, the tangential one
%                   % is needed for the return mapping algorithm.
%                   penalty_normal = IFpenalty_normal;
%                   penalty_tangent = IFpenalty_tangential;
%                 otherwise
%                   error('MATLAB:XFEM','Unvalid Sliding-ID');
%               end;
%               
%               % get residual contributions vectors, but choose betweeb the
%               % one- and two-integral formulation
%               switch IFintegral
%                 case 1  % one integral
%                   [res_penalty_normal res_penalty_tangent id tgappl ttrac f_trial] = ...
%                     get_ele_residual_penalty_alternative(xcoords,ycoords, ...
%                     seg_cut_info(i,e),INTERFACE_MAP(i).endpoints, ...
%                     id_dof(elenodes,:),id_eqns(elenodes,:),totaldis, ...
%                     penalty_normal,penalty_tangent,IFyieldstress, ...
%                     deltaload,IFsliding_switch);
%                 case 2  % two integrals
%                   [res_penalty_normal res_penalty_tangent id] = ...
%                     get_ele_residual_penalty(xcoords,ycoords, ...
%                     seg_cut_info(i,e),INTERFACE_MAP(i).endpoints, ...
%                     id_dof(elenodes,:),id_eqns(elenodes,:),totaldis);
%                 otherwise
%                   error('MATLAB:XFEM:UnvalidID', ...
%                     'Unvalid number of integrals');
%               end;
% 
%               % store plastic contribution to tangential gap
%               seg_cut_info(i,e).tgappl = tgappl;
%               
%               % store scalar values of current tangential traction at each
%               % gauss point
%               seg_cut_info(i,e).ttrac = ttrac;
%               
%               % store evaluated flow rules of trial state at both gauss
%               % points
%               seg_cut_info(i,e).f_trial = f_trial;
%               
%               % compute ele_residual_constraint
%               ele_residual_constraint = res_penalty_normal ...  % normal traction
%                                       + res_penalty_tangent;    % tangential traction
%                                                                  
%               % assemble into global residual
%               for a = 1:length(id)
%                 if id(a) ~= 0
%                   residual(id(a)) = residual(id(a)) + ele_residual_constraint(a);
%                 end;
%               end;
%             case 2  % Nitsche's mehtod
%             otherwise
%               error('MATLAB:XFEM','Unvalid Method-ID');
%           end;
%         end;
%       end;
%     end;
% 
%     residual = residual - big_force_loadstep;
%     
%     % clear some temporary variables
%     clear eleID elenodes xcoords ycoords res_constant_tangent i e a id ...
%       ele_residual_constraint constant_tangent res_penalty_normal ...
%       penalty_normal res_penalty_tangent penalty_tangent ...
%       normaltraction tangentialtraction;
    % ------------------------------------------------------------------- %
    %% UPDATE DISPLACEMENTS
%     % define a line serch parameter 'ls_par'
%     ls_par = 1;
%     
%     % check, if a line search is necessary
%     stol = 0.5;   % The value 'stol = 0.5' is recommended in 'Matthies1979'
%     G = totaldis * residual;
%     G0 = totaldis * residual_old;
%     if abs(G) > stol * abs(G0)
%       % compute the line search parameter 'ls_par'
%       ls_par = getlsparameter(stol,totaldis * residual_old, ...
%         totaldis * residual,totaldis);
%     end;
% 
%     deltanewton = ls_par * deltanewton;
    
    deltaload = deltaload + deltanewton;
    totaldis = totaldis + deltanewton;
    % ------------------------------------------------------------------- %
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
      ndisp(nnode,doff) = totaldis(i);
    end;
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
    end;
    
    % clear some temporary variables
    clear i j nnode doff grain;
    % ------------------------------------------------------------------- %
    %% COMPUTE SOME NORMS TO PREPARE CONVERGENCE CHECK      
    % store first displacement increment in order to use a relative
    % increment measure for the convergence check. The norms are normalized
    % to these values in each Newton-step.
    if iter  == 1
      deltanewton_null = norm(deltanewton);
      res_null = norm(residual);
      energ_null = residual' * deltanewton;
    end;

    % convergence check
    % The newton-scheme is assumed as converged, when the relative norms of
    % the residual, the displacement-increment and the energy-increment are
    % smaller than a given convergence tolerance 'IFconvtol' and there are 
    % no changes in the sliding states in comparison to the previous 
    % iteration step.
    
    % compute normalized norm of residual
    res_norm = norm(residual) / res_null;

    % compute normalized norm of diplacement inrement
    deltanewton_norm = norm(deltanewton) / deltanewton_null;
    
    % compute normalized norm of energie increment
    energ_norm = residual' * deltanewton / energ_null;

    % assume, that slidestates are converged
    slidestateconv = 'true';   % assume, that slidestates are converged

    % check, if slidestates are converged (only for nonlinear cases)
    if (IFsliding_switch == 2) || (IFsliding_switch == 3)    % only for perfect plasticity or frictional contact
      % loop over all interfaces 'i'
      for i=1:size(seg_cut_info,1)
        % loop over all elements
        for e=1:size(seg_cut_info,2)
          % only cut elements
          if seg_cut_info(i,e).elemno ~= -1
            if seg_cut_info(i,e).slidestateconv == 0
              % set 'slidestateconv'
              slidestateconv = 'false';
              
              % make sure, that the for-loops will be aborted as soon as
              % possible
              i=size(seg_cut_info,1);
              break;
            end;
          end;
        end;
      end;
    end;
    
%     % assume, that contact states are converged
%     contactstateconv = 'true';
% 
%     % check, if contact states are converged
%     if IFsliding_switch == 4  % only for case with frictionless contact
%       % loop over all interfaces 'i'
%       for i=1:size(seg_cut_info,1)
%         % loop over all elements
%         for e=1:size(seg_cut_info,2)
%           % only cut elements
%           if seg_cut_info(i,e).elemno ~= -1
%             if seg_cut_info(i,e).contactstateconv == 0
%               % set 'slidestateconv'
%               contactstateconv = 'false';
%               
%               % make sure, that the for-loops will be aborted as soon as
%               % possible
%               i=size(seg_cut_info,1);
%               break;
%             end;
%           end;
%         end;
%       end;
%     end;
    
    % plot evolution of slip-stick-area during Newton steps
    % Uncomment the following call to show possible oszillations during the
    % Newton iterations.
%     if timestep == 28
%     plotslidestateevolutionNewton;  % call a subroutine
%     end

    % build a formatted string with some information about the current
    % newton step
    iterdata = sprintf('  Iter %d  Res: %.4e  Delta: %.4e  Energ: %.4e', ...
      iter,res_norm,deltanewton_norm,energ_norm);
%     iterdata = sprintf('  Iter %d  Res: %.4e  Delta: %.4e  Energ: %.4e  Slide: %s \t Cont: %s' , ...
%       iter,res_norm,deltanewton_norm,energ_norm,slidestateconv,contactstateconv);
    

    % print some information about the current newton step
    disp(iterdata);
    
    % clear some temporary variables
    clear i j iterdata;
    % ------------------------------------------------------------------- %
    %% CHECK CONVERGENCE OF NEWTON-RAPHSON-SCHEME 
    if (res_norm < IFconvtol) && (deltanewton_norm < IFconvtol) && ...
      (energ_norm < IFconvtol);% && (strcmp(slidestateconv,'true') == 1) && ...
%       (strcmp(contactstateconv,'true') == 1)
      break;  % exit iteration loop (while-loop), if convergence is achieved
    end;

    % check, if number of maximum iterations is reached
    if iter >= IFmaxiter
      error('MATLAB:XFEM:main_xfem',...
          'Newton did not converge in %d iterations.', IFmaxiter);
    end;

    % increase the iteration index 'iter'
    iter = iter + 1;
    % ----------------------------------------------------------------- %
  end;
  
  % clear some temporary variables
  clear deltanewton deltaload iter delta_norm res_norm energ_norm ...
    deltanewton_null res_null energ_null slidestateconv;
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
    ndisp(nnode,doff) = totaldis(i);
  end;
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
  end;
  
  % clear some temporary variables
  clear i j nnode doff grain;
  % --------------------------------------------------------------------- %
  %% UPDATE AT END OF EACH LOAD-STEP  
% clear x_def y_def;
% x_def = [];
% y_def = [];
% 
% % get coordinates of deformed mesh
% for i = 1:numnod
%     x_def(i) = x(i) + dis(2*i-1);
%     y_def(i) = y(i) + dis(2*i); 
% end
  
  % updates at the end of each load step, that are specific with respect to
  % the sliding case
  switch IFsliding_switch
    case 0  % fully tied case
      % no updates necessary
    case 1  % frictionless sliding
      % no updates necessary
      
%       figure(51)
%       hold on;
%       plot(x_def(881)-x(881),x_def(891)-x(891),'*') ;
%     %   axis([-0.012 0.008 -0.0005 0.0025]);
%       hold off;

%       figure(52)
%       hold on;
%       plot(x_def(1),y_def(1),'*r'); 
%       hold off;
      
      % figure(51)
%   hold on;
%   plot(x_def(3361)-x(3361),x_def(3381)-x(3381),'*') 
% %   axis([-0.012 0.008 -0.0005 0.0025]);
%   hold off;
%   
%   figure(52)
%   hold on;
%   plot(x_def(1),y_def(1),'*') 
%   hold off;
    case 2  % perfect plasticity
      %% update plastic contribution to plastic gap and tangential tractions
      % loop over interfaces 'i'
      for i=1:size(seg_cut_info,1)
        % loop over cut elements 'e'
        for e=1:size(seg_cut_info,2)
          % only, if 'e' is cut by 'i'
          if seg_cut_info(i,e).elemno ~= -1
            % current (converged) state is the converged step for the next
            % load step

            % plastic contribution to tangential gap
            seg_cut_info(i,e).tgapplconv = seg_cut_info(i,e).tgappl;

            % reset current gap
            seg_cut_info(i,e).tgappl = [0 0];

            % current (converged) tangential traction is the converged step 
            % for the next load step
            seg_cut_info(i,e).ttracconv = seg_cut_info(i,e).ttrac;

            % reset current tangential traction at the two gauss points
            seg_cut_info(i,e).ttrac = [0 0];
          end;
        end;
      end;
      % ----------------------------------------------------------------- %
      %% visualize 'slidestate'-flags
%{
      % create a new figure
      if timestep == 1
        figure(30)
      end;
      set(0,'CurrentFigure',30)
      hold on;

      % plot interfaces
      % loop over all interfaces
      for i=1:size(seg_cut_info,1)
        % loop over all elements
        for e=1:size(seg_cut_info,2)
          % only cut elements
          if seg_cut_info(i,e).elemno ~= -1
            % get 2 points, that determine the subsegment
            if all(size(seg_cut_info(i,e).xint) == [2 2])       % no triple junction in 'e'
                p1 = seg_cut_info(i,e).xint(1,:);
                p2 = seg_cut_info(i,e).xint(2,:);
            elseif all(size(seg_cut_info(i,e).xint) == [1 2])   % triple junction in 'e'
                p1 = seg_cut_info(i,e).xint(1,:);

                % Get coordinates of nodes of element
                xep=zeros(1,3);
                yep=zeros(1,3);
                for m=1:3
                    jep = node(m,seg_cut_info(i,e).elemno); 
                    xep(m) = x(jep); 
                    yep(m) = y(jep);
                end

                % Second endpoint of segment is also end point of
                % interface --> check, which one of the two endpoints
                % of the interface lies in element 'e'

                % get first endpoint
                endpoint = INTERFACE_MAP(i).endpoints(1,:); 

                % check, if it is in the element 'e'
                inside = polygon_contains_point_2d ( 3, [xep;yep], endpoint );

                if inside       % endpoint is in element
                    p2 = endpoint;
                else            % endpoint is not in element
                    p2 = INTERFACE_MAP(i).endpoints(2,:);
                end;
            end;
            % now, p1 and p2 are the two nodes, that determine the
            % subsegment
            xcoord = [p1(1) p2(1)];
            ycoord = [p1(2) p2(2)];

            % set color depending on current slidestate
            stylecell = {'b','r'};

            if any(seg_cut_info(i,e).f_trial > 0)
              index = 2;
            else
              index = 1;
            end;
            
%             % get 'index' for stylecell, depending on the 'slidestate'
%             index = seg_cut_info(i,e).slidestate + 1;

            % get current normalized pseudotime of entire simulation
            timecoord = (timestep) / (length(time) - 1);
            
            % check, if it is a horizontal or a vertical interface
            if abs(xcoord(1) - xcoord(2)) < abs(ycoord(1)-ycoord(2))
              % vertical interface
              plot([timecoord timecoord],ycoord,stylecell{index},'LineWidth',3.0);   
              xlabel('normalized pseudo-time');
              ylabel('y-coordinate');
            else
              % horizontal interface
              plot(xcoord,[timecoord timecoord],stylecell{index},'LineWidth',3.0);       
              xlabel('x-coordinate');
              ylabel('normalized pseudo-time');
            end;
          end;
        end;
      end;

      % edit figure
      % legend('stick','slip');
      title('blue = stick, red = slip');
      hold off;
%}
      % ----------------------------------------------------------------------- %
      
%       % plot
%   figure(50)
%   hold on;
% %   dofD = id_eqns(1,1)
% %   dofC = id_eqns(881,1)
%   plot(x(881),y(881),'*g');
%   plot(x(891),y(891),'*r');
%   hold off;

% x_def = [];
% y_def = [];
% 
% % get coordinates of deformed mesh
% for i = 1:numnod
%     x_def(i) = x(i) + dis(2*i-1);
%     y_def(i) = y(i) + dis(2*i); 
% end
%   
%       figure(51)
%       hold on;
%       plot(x_def(881)-x(881),x_def(891)-x(891),'*') ;
%     %   axis([-0.012 0.008 -0.0005 0.0025]);
%       hold off;
% 
%       figure(52)
%       hold on;
%       plot(x_def(881),y_def(881),'*r'); 
%       hold off;

% figure(51)
%   hold on;
%   plot(x_def(3361)-x(3361),x_def(3381)-x(3381),'*') 
% %   axis([-0.012 0.008 -0.0005 0.0025]);
%   hold off;
%   
%   figure(52)
%   hold on;
%   plot(x_def(1),y_def(1),'*') 
% %   axis([-0.012 0.008 -0.0005 0.0025]);
%   hold off;
    case 3  % frictional contact (Coulomb)
      warning('MATLAB:XFEM:main_xfem',...
        'There exists no code for frictional contact (Coulomb), yet.')
    case 4  % frictionless contact (only opening contact)
      % no updates necessary
    otherwise
      warning('MATLAB:XFEM:main_xfem',...
        'Unvalid slidingID. Choose valid ID or add additional case to switch-case-structure')
  end;  
  
  % clear some temporary variables
  clear x_def y_def i;
  % --------------------------------------------------------------------- %
  %% LOAD STEPPING LOOP (END)
end;
% ----------------------------------------------------------------------- %
%% RE-ASSEMBLE GLOBAL DISPLACEMENT VECTOR
% Reassemble displacement vector  - Enriched nodes need to have their
% degrees of freedom added to end up representing a total displacement.
% The extra degrees of freedom should only be added if the node is enriched
% with the grain in which it resides.

% This re-assembly process strictly uses the formulation in 'Simone2006'.

% clear some variables, which might be defined due to Nitsche's method
clear ndisp x_def y_def;

% ndisp(i,:) are all of the solutions (base and enriched) at node i
ndisp = zeros(numnod,6);

% dis is a vector with traditional FEM numbering, and the final solution
% for each node. It stores the nodal displacements
dis = zeros(2*numnod,1);

for i = 1:numeqns
  [nnode,doff] = find(id_eqns == i);
  ndisp(nnode,doff) = totaldis(i);
end

% copy 'ndisp' before modifying it, since the original one is needed for
% postprocessing the stresses
orig_ndisp = ndisp;

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

% clear some temporary variables
clear i j nnode doff grain;
% ----------------------------------------------------------------------- %
%% POST PROCESS: STRESSES
disp('postprocessing: stresses ...');


% compute stresses 'stress' and strains 'strain' at center of each element
% Structure of 'stress':
%   Dimension: numelex6xnumgrains
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


% von-Mises-stresses at center of each element
% Structure of 'stressvonmises':
%   Dimension: numelex4xnumgrains
%   Index 1:    'stress' for element 'e'
%   Index 2:    column 1    global element ID
%               column 2    x-coordinate of element centroid
%               column 3    y-coordinate of element centroid
%               column 4    von-Mises-stress at element centroid
%   Index 3:    ID of grain, to which these values belong to
stressvonmises = zeros(numele,4,maxngrains);


% loop over all elements
for e=1:numele
%     [stresse] = post_process(node,x,y,e,dis);
%     stress(e,1:6) = stresse;
  
  % compute stress and strain in current element
  [straine,stresse] = post_process_better(node,x,y,e,dis,orig_ndisp, ...
    id_dof,cutlist,maxngrains,INT_INTERFACE);
  
  % assign 'stresse' to 'stress'
  stress(e,1:6,:) = stresse;
  
  % assign 'straine' to 'strain'
  strain(e,1:6,:) = straine;
  
  % compute von-Mises-stress
  stressvonmises(e,1:3,:) = stresse(1,1:3,:);
  for i=1:maxngrains
    stressvonmises(e,4,i) = sqrt((stresse(1,4,i))^2 + (stresse(1,5,i))^2 - ...
      stresse(1,4,i) * stresse(1,5,i) + 3 * stresse(1,6,i)^2);
  end;
end


% clear some temporary variables
clear stresse straine maxstress_vec minstress_vec i e f j;
% ----------------------------------------------------------------------- %
%% POST PROCESS: TRACTIONS AT INTERFACE
disp('postprocessing: tractions ...');
  %% direct evaluation
  % computing the internal forces in the interface depends on the method of
  % constraint enforcing
  switch IFmethod 
    case 0  % Lagrange multipliers
      % extract vector with lagrange multipliers from 'totaldis'
      lagmult = totaldis(old_size+1:old_size + 2*multipliers);

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

      for i = 1:size(seg_cut_info,1)     % for every interface
        for e = 1:size(seg_cut_info,2) % for every cut element in that interface
          % initialize variable for lagrange multiplier in 'seg_cut_info'
          seg_cut_info(i,e).lagmult = [];

          if seg_cut_info(i,e).elemno ~= -1
            % get some element data
            eleID = seg_cut_info(i,e).elemno;
            elenodes = node(:,eleID);
            xcoords = x(elenodes);
            ycoords = y(elenodes);
            
            % distinguish between sliding cases
            switch IFsliding_switch
              case 0  % fully tied case
                penalty_normal = IFpenalty_normal;
                penalty_tangent = IFpenalty_tangential;
              case 1  % frictionless sliding
                penalty_normal = IFpenalty_normal;
                penalty_tangent = 0;
              case 2  % perfect plasticity
              otherwise
                error('MATLAB:XFEM:UnvalidID', ...
                  'Unvalid sliding ID');
            end;
            
            % choose between one- and two-integral formulation
            switch IFintegral
              case 1  % one integral
                [normaltraction tangentialtraction] = ...
                  get_lag_mults_for_penalty_alternative(xcoords, ...
                    ycoords,seg_cut_info(i,e), ...
                    INTERFACE_MAP(i).endpoints,id_dof(elenodes,:), ...
                    id_eqns(elenodes,:),totaldis);
              case 2  % two integrals
                [normaltraction tangentialtraction] = ...
                  get_lag_mults_for_penalty(xcoords, ...
                    ycoords,seg_cut_info(i,e), ...
                    INTERFACE_MAP(i).endpoints,id_dof(elenodes,:), ...
                    id_eqns(elenodes,:),totaldis');
              otherwise
                error('MATLAB:XFEM:UnvalidID', ...
                  'Unvalid number of integrals');
            end;
            
            normaltraction = penalty_normal * normaltraction;
            tangentialtraction = penalty_tangent * tangentialtraction;
            
            seg_cut_info(i,e).lagmult = [normaltraction tangentialtraction];

            lagmult = [lagmult normaltraction tangentialtraction];
          end;
        end;
      end;

      % clear some temporary variables
      clear pn_nodes neg_g pos_g i e parent_el;
    case 2  % Nitsche's method
      % Compute Lagrange multipliers via 
      %                       'lambda = alpha * [[u]] - <sigma>*normal'

      for i = 1:size(seg_cut_info,1)     % for every interface
        for e = 1:size(seg_cut_info,2) % for every cut element in that interface
          % initialize variable for lagrange multiplier in 'seg_cut_info'
          seg_cut_info(i,e).lagmult = [];

          if seg_cut_info(i,e).elemno ~= -1

            % get stabilization parameter
            penalty = seg_cut_info(i,e).nitsche;
            
            % get current element ID
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
              INTERFACE_MAP(i).endpoints,penalty,totaldis');

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
  
  % clear some temporary variables
  clear i j normaltraction tangentialtraction eleID elenodes xcoords ...
    ycoords penalty_normal penalty_tangent normaltraction ...
    tangentialtraction;
  % --------------------------------------------------------------------- %
  %% domain integral method 
%   % Since the constraint enforcement methods represent the tractions at the
%   % interface not very accurate, the domain integral method based on [Dolbow,
%   % J. and Harari, I.: An efficient finite element method for embedded
%   % interface problems. Int. J. Numer. Meth. Engng. 2009 (78): 229-252] can
%   % be used to get a better approximation of the flux across the interface.
%   %
%   % Actually, almost all quantities, that are needed for the following
%   % computation, are already known, since they are part of 'bigk' or
%   % 'big_force'. The right entries have to be extracted from those
%   % structures.
% 
%   % loop over all interfaces 'i'
%   for i=1:size(seg_cut_info,1)
%     % loop over all elements 'e' of that interface
%     for e=1:size(seg_cut_info,2)
%       % only cut elements
%       if seg_cut_info(i,e).elemno ~= -1
%         % get global element ID of current element
%         eleID = seg_cut_info(i,e).elemno;
% 
%         % get nodes of current element
%         elenodes = node(:,eleID);
% 
%         % get DOFs of current element nodes 'elenodes'
%         DOFs = id_eqns(elenodes,:);
% 
%         % get enriching grains for the element's nodes
%         enr_grains = id_dof(elenodes,:);
% 
%         % extract force values from 'big_force'
%   %       domint_force(1) = big_force(DOFs(1,1));
%   %       domint_force(2) = big_force(DOFs(1,2));
%   %       domint_force(3) = big_force(DOFs(2,1));
%   %       domint_force(4) = big_force(DOFs(2,2));
%   %       domint_force(5) = big_force(DOFs(3,1));
%   %       domint_force(6) = big_force(DOFs(3,3));
%   %       domint_dis(1) = totaldis(DOFs(1,1));
%   %       domint_dis(2) = totaldis(DOFs(1,2));
%   %       domint_dis(3) = totaldis(DOFs(2,1));
%   %       domint_dis(4) = totaldis(DOFs(2,2));
%   %       domint_dis(5) = totaldis(DOFs(3,1));
%   %       domint_dis(6) = totaldis(DOFs(3,2));
%   %       
%   %       if DOFs(1,3) ~= 0 && DOFs(1,4) ~= 0
%   %         domint_force(7) = big_force(DOFs(1,3));
%   %         domint_force(8) = big_force(DOFs(1,4));
%   %         domint_dis(7) = totaldis(DOFs(1,3));
%   %         domint_dis(8) = totaldis(DOFs(1,4));
%   %       end;
%   %       
%   %       if DOFs(2,3) ~= 0 && DOFs(2,4) ~= 0
%   %         domint_force(9) = big_force(DOFs(2,3));
%   %         domint_force(10) = big_force(DOFs(2,4));
%   %         domint_dis(9) = totaldis(DOFs(2,3));
%   %         domint_dis(10) = totaldis(DOFs(2,4));
%   %       end;
%   %       
%   %       if DOFs(3,3) ~= 0 && DOFs(3,4) ~= 0
%   %         domint_force(11) = big_force(DOFs(3,3));
%   %         domint_force(12) = big_force(DOFs(3,4));
%   %         domint_dis(11) = totaldis(DOFs(3,3));
%   %         domint_dis(12) = totaldis(DOFs(3,4));
%   %       end;
%   %       
%   %       if DOFs(1,5) ~= 0 && DOFs(1,6) ~= 0
%   %         domint_force(13) = big_force(DOFs(1,5));
%   %         domint_force(14) = big_force(DOFs(1,6));
%   %         domint_dis(13) = totaldis(DOFs(1,5));
%   %         domint_dis(14) = totaldis(DOFs(1,6));
%   %       end;
%   %       
%   %       if DOFs(2,5) ~= 0 && DOFs(2,6) ~= 0
%   %         domint_force(15) = big_force(DOFs(2,5));
%   %         domint_force(16) = big_force(DOFs(2,6));
%   %         domint_dis(15) = totaldis(DOFs(2,5));
%   %         domint_dis(16) = totaldis(DOFs(2,6));
%   %       end;
%   %       
%   %       if DOFs(3,5) ~= 0 && DOFs(3,6) ~= 0
%   %         domint_force(17) = big_force(DOFs(3,5));
%   %         domint_force(18) = big_force(DOFs(3,6));
%   %         domint_dis(17) = totaldis(DOFs(3,5));
%   %         domint_dis(18) = totaldis(DOFs(3,6));
%   %       end;
% 
% 
%         for j=1:size(DOFs,2)
%           if DOFs(1,j) ~= 0
%             domint_force(j) = big_force(DOFs(1,j));
%           end;
%           if DOFs(2,j) ~= 0
%             domint_force(6+j) = big_force(DOFs(2,j));
%           end;
%           if DOFs(3,j) ~= 0
%             domint_force(12+j) = big_force(DOFs(3,j));
%           end;
%         end;
% 
%         % rearrange 'DOFs' as a vector
%         DOFs_vec = [DOFs(1,:) DOFs(2,:) DOFs(3,:)];
% 
%         % extract stiffness values form 'bigk'
%         for j=1:length(DOFs_vec)
%           for k = 1:length(DOFs_vec)
%             if DOFs_vec(j) ~= 0 && DOFs_vec(k) ~= 0
%               domint_stiff(j,k) = bigk(DOFs_vec(j),DOFs_vec(k));
%             end;
%           end;
%         end;
% 
%         % extract displacement values from 'totaldis'
%         for j=1:size(DOFs,2)
%           if DOFs(1,j) ~= 0
%             domint_dis(j) = totaldis(DOFs(1,j));
%           end;
%           if DOFs(2,j) ~= 0
%             domint_dis(6+j) = totaldis(DOFs(2,j));
%           end;
%           if DOFs(3,j) ~= 0
%             domint_dis(12+j) = totaldis(DOFs(3,j));
%           end;
%         end;
% 
%         % call a subroutine to compute the traction in this element
%         seg_cut_info(i,e).domint = evaluatedomainintegral(elenodes,DOFs, ...
%           id_dof,seg_cut_info(i,e),x(elenodes),y(elenodes),domint_force, ...
%           domint_stiff,domint_dis);
%       end;
%     end;
%   end;
% 
%   % clear some temporary variables
%   clear i e eleID elenodes DOFs;
%   % --------------------------------------------------------------------- %
%% POST PROCESS: PREPARE showdeform2
x_def = [];
y_def = [];

% get coordinates of deformed mesh
for i = 1:numnod
    x_def(i) = x(i) + dis(2*i-1);
    y_def(i) = y(i) + dis(2*i); 
end

% clear some temporary variables
clear i;
% ----------------------------------------------------------------------- %
%% FINISH SOLVING PROCESS
%  disp('saving to results file ...');

%  save my_results_files.mat x y node dis totaldis numele numnod cutlist...
%      INT_INTERFACE SUBELEM_INFO SUBELEMENT_GRAIN_MAP id_eqns id_dof...
%      X Y CONN stress

disp('Solving and postprocessing done.');

