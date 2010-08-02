% Input File 'inp_patchtest_24_6.m'
%
% Here, you can define all parameters to configure the simulation.
%
%**************************************************************************
% GIVE A SHORT DESCRIPTION OF THE EXAMPLE
%**************************************************************************
% Rectangle. Length x Heigth = 16 x 4. Mesh 24 x 6. Same Material in every
% grain (E = 1000.0, nue = 0.3). Constant load on right side, pulling in
% x-direction. Left side fixed. 3 grains.
%**************************************************************************
%
% To set up a new example, build it in this file, so that all IDs are
% strored in this file. Save it to an own file, afterwards.
%

% Author: Matthias Mayr (04/2010)

%--------------------------------------------------------------------------
% PARAMETERS FOR BACKGROUND MESH
% Set some parameters to specifiy the mesh and the geometry
% 
% mesh structure: 'IFmeshstructure'
% ID    Description
% 0     structured
% 1     unstructured
% 2     read mesh from gmsh-mesh-file '*.msh'
IFmeshstructure = 0;
%
% Shape of geometry: 'IFshapegeometryID'
% ID    Description
% 0     rectangular
%
IFshapegeometryID = 0;
%
% Give length and heigth of rectangle
IFlength = 16;
IFheight = 4;
%
% Give number of line divisions in x- and y-direction
IFnldivx = 24;
IFnldivy = 6;
%
% filename of msh-file withput file extension '.msh'
% (if reading mesh from gmsh-msh-file)
IFfilename_msh_file = 'quarterring_gmsh';      % NO FILE EXTENSION '.msh'
%--------------------------------------------------------------------------
% PARAMETERS FOR INTERFACES
% Set some parameters to specify the interfaces (boundaries of the grains)
%
% Choose one of the datasets for p in 'comp_geo/vdata_multi.m'
%
IFdatasetp = 4;%4;%11;%19;
%--------------------------------------------------------------------------
% BOUNDARY CONDITIONS
% Dirichlet Boundary Conditions (DBCs) and Neumann Boundary Conditions
% (NBCs) are applied separately. So, you have to chose a set of BCs for
% both.
%
% Dirichlet BCs
% ID    Description
% 0     bc_conv4_DBC.m
% 1     bc_conv1_DBC.m
% 2     bc_conv7_DBC.m
% 3     bc_conv9_DBC.m
% 4     bc_conv10_DBC.m
% 5     bc_conv11_DBC.m
% 6     multi1_DBC.m
% 7     frictionless_sliding1_24_6_DBC.m
% 8     frictionless_sliding1_72_18_DBC.m
% 9     quartercircle_gmsh_DBC.m
% 10    patchtest_8_2_DBC.m
% 11    patchtest_24_6_DBC.m
IFDirichletBCs = 11;
%
% Neumann BCs
% ID    Filename            Description
% 0     bc_conv4_NBC.m
% 1     bc_conv1_NBC.m
% 2     bc_conv7_NBC.m
% 3     bc_conv9_NBC.m
% 4     bc_conv10_NBC.m
% 5     bc_conv11_NBC.m
% 6     multi1_NMC.m
% 7     frictionless_sliding1_24_6_NBC.m
% 8     frictionless_sliding1_72_18_NBC.m
% 9     quartercircle_gmsh_NBC.m
% 10    patchtest_8_2_NBC.m
% 11    patchtest_24_6_NBC.m
IFNeumannBCs = 11;
%
% method of giving NBCs
% ID    Description
% 0     nodal forces (integration done by user, only not-enriched nodes)
% 1     tractions given as functions
IFneumann = 0;
%--------------------------------------------------------------------------
% MATERIAL PROPERTIES
% Set an ID 'IFMatSet' to chose a set of material properties from material 
% database 'preprocess\MaterialProperties.m'
% ID    Description
% 0     24 grains with same material properties (nue = 0.3, E = 1000.0)
% 1     Two grains (nue1 = 0.0, nue2 = 0.3, E1 = E2 = 1000.0)
% 2     24 grains with different material properties
% 3     24 grains with same material properties (nue = 0.0, E = 1000.0)
IFMatSet = 3;
%--------------------------------------------------------------------------
% METHOD OF ENFORCING CONSTRAINTS AT THE INTERFACE
% Set an ID to choose the method, by which the constrains shall be enforced
% at the interface
% ID    Method
% 0     Lagrange Multipliers (piecewise constant)
% 1     Penalty-Method
% 2     Nitsche's Method
IFmethod = 1;
%
% Set Penalty-Parameter
IFpenalty_normal      = 5e+10;
IFpenalty_tangential  = 5e+10;
%
% Nitsche Parameter
IFnitsche = 1.0e+3;
%
% Choose a penalty variant: One or two integrals
% ID    Number of integrals
% 1     One integral (alpha ~1/h)
% 2     Two integrals (alpha ~1/h^2)
IFintegral = 1;
%--------------------------------------------------------------------------
% SLIDING PARAMETERS
% Set an ID to indicate, how sliding should be treaten: 'IFsliding_switch'
% ID    Description
% 0     no sliding at all (full constraint)
% 1     frictionless sliding
% 2     perfect plasticity with shear yield stress
% 3     frictional sliding with Coulomb's friction
%
IFsliding_switch = 2; 
% 
% Set a yield stress for plasticity
IFyieldstress = 1.3;
%--------------------------------------------------------------------------
% SOLVER PREFERENCES
% You can choose between an explicit solver and an implicit solver via a
% Newton-Raphson-scheme. For the implicit one, you have to set some
% additional preferences
%
% Solver Type
% ID    Description
% 0     explicit
% 1     implicit (Newton-Raphson-scheme)
IFSolverType = 1;
%
% Maximum number of iterations 'IFmaxiter' (only for implicit solver)
IFmaxiter = 50;
%
% convergence criteria: increment of displacement < 'IFconvtol' ???
IFconvtol = 1.0e-7;
%
% vector with pseudo-time-steps (always between '0' and '1')
IFtime = linspace(0,1,501);  % vector creation without 'linspace'-command
                           % possible, but first element has to be '0'
%--------------------------------------------------------------------------
% THE PARAMETER LIST ENDS HERE. DO NOT TOUCH ANY CODE BEYOND THIS LINE !!!
%--------------------------------------------------------------------------

% print parameters in console
disp(['IFmeshstructure:         ' num2str(IFmeshstructure)]);
disp(['IFshapegeometryID:       ' num2str(IFshapegeometryID)]);
disp(['IFlength:                ' num2str(IFlength)]);
disp(['IFheight:                ' num2str(IFheight)]);
disp(['IFnldivx:                ' num2str(IFnldivx)]);
disp(['IFnldivy:                ' num2str(IFnldivy)]);
disp(['IFdatasetp:              ' num2str(IFdatasetp)]);
disp(['IFDirichletBCs:          ' num2str(IFDirichletBCs)]);
disp(['IFNeumannBCs:            ' num2str(IFNeumannBCs)]);
disp(['IFMatSet:                ' num2str(IFMatSet)]);
disp(['IFsliding_switch:        ' num2str(IFsliding_switch)]);
disp(['IFmethod:                ' num2str(IFmethod)]);
disp(['IFpenalty_normal:        ' num2str(IFpenalty_normal)]);
disp(['IFpenalty_tangential:    ' num2str(IFpenalty_tangential)]);
disp(['IFnitsche:               ' num2str(IFnitsche)]);
disp(['IFSolverType:            ' num2str(IFSolverType)]);
disp(['IFmaxiter:               ' num2str(IFmaxiter)]);
disp(['IFconvtol:               ' num2str(IFconvtol)]);

%--------------------------------------------------------------------------
% create filenames for input files from type '*.mat'
filename1 = fullfile(pwd, 'comp_geo', 'xfeminputdata_comp_geo.mat');
filename2 = fullfile(pwd, 'preprocess', 'xfeminputdata_preprocess.mat');
filename3 = fullfile(pwd, 'XFEM', 'xfeminputdata_xfem.mat');

% save input parameters to mat-files into specific directories
save(filename1, 'IFmeshstructure', 'IFshapegeometryID', 'IFlength', ...
    'IFheight', 'IFnldivx', 'IFnldivy', 'IFdatasetp');  % for 'comp_geo'
save(filename2, 'IFDirichletBCs', 'IFNeumannBCs', 'IFMatSet'); % for 'preprocess'
save(filename3, 'IFsliding_switch','IFmethod','IFpenalty_normal', ...
  'IFpenalty_tangential','IFnitsche','IFSolverType','IFmaxiter', ...
  'IFconvtol');                    % for 'XFEM'

% clear workspace
clear IFmeshstructure IFshapegeometryID IFlength IFheight IFnldivx ...
    IFnldivy IFdatasetp IFDirichletBCs IFNeumannBCs IFMatSet ...
    IFsliding_switch IFmethod IFpenalty IFnitsche IFSolverType ...
    IFmaxiter IFconvtol;
clear filename1 filename2 filename3;