% Preprocessing for the Lagrange Multiplier XFEM code.  Given a mesh
% geometry, already divided into sub-meshes via the comp_geo routines, this
% routine defines boundary conditions and material properties.
% Jessica Sanders, Duke University, 2010

% ----------------------------------------------------------------------- %

% LOAD MODEL DATA

load cstruct4.mat

% Get the number of total subelements
num_sub_elems = size(CONN,2);

% ----------------------------------------------------------------------- %
   
% BOUNDARY CONDITIONS  - FORCE AND DISPLACEMENTS

f = 0;  % Some boundary conditions require the parameter, f
[force, dispbc, ubar, num_enr_surf, enr_surf, bc_enr]...
    = applybcs(x,y,numnod,beam_l,beam_h,f);

% Do the enriched boundaries have sliding? 1 = yes, 0 = no
sliding_switch = 0;


% ----------------------------------------------------------------------- %

% GRAIN MATERIAL INFORMATION

global GRAININFO_ARR
GRAININFO_ARR = struct('grain_no', 0, 'num_elems',0,'poisson',0,'youngs',0);

% assign material properties to each grain
poissons = [0.3 0.3 0.3 0.3];
youngs = [1000.0 1000.0 1000.0 1000.0];
neg = sum(elemgrainmap);    % number of elements in each grain

for i = 1:maxngrains
    GRAININFO_ARR(i).grain_no = i;
    GRAININFO_ARR(i).num_elems = neg(i);
    GRAININFO_ARR(i).poisson = poissons(i);
    GRAININFO_ARR(i).youngs = youngs(i);
end

%Set a tolerance size (best to base it on the mesh size)
tol = 0.00001;

% Save all variables into a input file

save beam_bending_example4.mat