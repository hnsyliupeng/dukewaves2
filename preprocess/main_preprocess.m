% Preprocessing for the Lagrange Multiplier XFEM code.  Given a mesh
% geometry, already divided into sub-meshes via the comp_geo routines, this
% routine defines boundary conditions and material properties.
% Jessica Sanders, Duke University, 2010

% ----------------------------------------------------------------------- %

% Load input parameters
load xfeminputdata_preprocess.mat

% LOAD MODEL DATA
load my_new_mesh.mat

% Get the number of total subelements
num_sub_elems = size(CONN,2);

% ----------------------------------------------------------------------- %
   
% BOUNDARY CONDITIONS  - FORCE AND DISPLACEMENTS

% Applying of BCs is managed by 'applybcs.m'
f = 0;  % Some boundary conditions require the parameter, f
[force, dispbc, ubar, num_enr_surf, enr_surf, bc_enr,nodeNBC]...
    = applybcs(x,y,numnod,beam_l,beam_h,f);     % call 'applybcs.m'
    % force  ...    Neumann BCs
    % dispbc ...    knows, in which DOFs there is a DBC
    % ubar   ...    Gives the values of the DBCs

% ----------------------------------------------------------------------- %

% GRAIN MATERIAL INFORMATION

global GRAININFO_ARR
GRAININFO_ARR = struct('grain_no', 0, 'num_elems',0,'poisson',0,'youngs',0);

% assign material properties to each grain
MatProp = MaterialProperties(IFMatSet);     % Call Material Database
poissons = MatProp.poissons;
youngs = MatProp.youngs;
neg = sum(elemgrainmap);    % number of elements in each grain

for i = 1:maxngrains
    GRAININFO_ARR(i).grain_no = i;
    GRAININFO_ARR(i).num_elems = neg(i);
    GRAININFO_ARR(i).poisson = poissons(i);
    GRAININFO_ARR(i).youngs = youngs(i);
end

%Set a tolerance size (best to base it on the mesh size)
tol = 0.00001;

% clear temporary variables from workspace
clear MatProp;

% Save all variables into a input file
save my_new_mesh_with_BCs.mat

% print message into console
disp('Preprocessing finished. Variables saved to file "my_new_mesh_with_BCs.mat"');
