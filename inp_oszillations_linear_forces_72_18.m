% Input File 'inp_oszillations_linear_forces_72_18.m'
%
% Here, you can define all parameters to configure the simulation.
%
%**************************************************************************
% GIVE A SHORT DESCRIPTION OF THE EXAMPLE
%**************************************************************************
% Rectangle. Length x Heigth = 16 x 4. Mesh 4 x 1. Same Material in every
% grain (E = 1000.0, nue = 0.3). Hatshaped forces, pulling in
% x-direction. Left side fixed. Should show oszillations for Lagrange
% multipliers.
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
IFnldivx = 72;
IFnldivy = 18;
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
IFdatasetp = 13;
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
% 12    patchtest_72_18_DBC.m
% 13    patchtest_4_1_DBC.m
% 14    patchtest_392_18_DBC.m
IFDirichletBCs = 12;
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
% 12    patchtest_72_18_NBC.m
% 13    patchtest_4_1_NBC.m
% 14    patchtest_392_98_NBC.m
% 15    hatshaped_x_forces_72_18_NBC.m
% 16    parabolic_x_forces_72_18_NBC.m
IFNeumannBCs = 15;
%--------------------------------------------------------------------------
% MATERIAL PROPERTIES
% Set an ID 'IFMatSet' to chose a set of material properties from material 
% database 'preprocess\MaterialProperties.m'
% ID    Description
% 0     24 grains with same material properties (nue = 0.3, E = 1000.0)
% 1     Two grains (nue1 = 0.0, nue2 = 0.3, E1 = E2 = 1000.0)
% 2     24 grains with different material properties
% 3     24 grains with same material properties (nue = 0.0, E = 1000.0)
% 4     3 grains with different material properties
IFMatSet = 3;
%--------------------------------------------------------------------------
% METHOD OF ENFORCING CONSTRAINTS AT THE INTERFACE
% Set an ID to choose the method, by which the constrains shall be enforced
% at the interface
% ID    Method
% 0     Lagrange Multipliers (piecewise constant)
% 1     Penalty-Method
% 2     Nitsche's Method
IFmethod = 0;
%
% Set Penalty-Parameter
IFpenalty = 5.0e+5;
%
% Nitsche Parameter
IFnitsche = 0;
%--------------------------------------------------------------------------
% SLIDING PARAMETERS
% Set an ID to indicate, how sliding should be treaten: 'IFsliding_switch'
% ID    Description
% 0     no slidung at all (full constraint)
% 1     frictionless sliding
% 2     perfect plasticity with shear yield stress
% 3     frictional sliding with Coulomb's friction
%
IFsliding_switch = 0; 
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
IFSolverType = 0;
%
% Maximum number of iterations 'IFmaxiter' (only for implicit solver)
IFmaxiter = 25;
%
% convergence criteria: increment of displacement < 'IFconvtol' ???
IFconvtol = 1.0e-8;
%--------------------------------------------------------------------------
% THE PARAMETER LIST ENDS HERE. DO NOT TOUCH ANY CODE BEYOND THIS LINE !!!
%--------------------------------------------------------------------------

% print parameters in console
disp(['IFmeshstructure:     ' num2str(IFmeshstructure)]);
disp(['IFshapegeometryID:   ' num2str(IFshapegeometryID)]);
disp(['IFlength:            ' num2str(IFlength)]);
disp(['IFheight:            ' num2str(IFheight)]);
disp(['IFnldivx:            ' num2str(IFnldivx)]);
disp(['IFnldivy:            ' num2str(IFnldivy)]);
disp(['IFdatasetp:          ' num2str(IFdatasetp)]);
disp(['IFDirichletBCs:      ' num2str(IFDirichletBCs)]);
disp(['IFNeumannBCs:        ' num2str(IFNeumannBCs)]);
disp(['IFMatSet:            ' num2str(IFMatSet)]);
disp(['IFsliding_switch:    ' num2str(IFsliding_switch)]);
disp(['IFmethod:            ' num2str(IFmethod)]);
disp(['IFpenalty:           ' num2str(IFpenalty)]);
disp(['IFnitsche:           ' num2str(IFnitsche)]);
disp(['IFSolverType:        ' num2str(IFSolverType)]);
disp(['IFmaxiter:           ' num2str(IFmaxiter)]);
disp(['IFconvtol:           ' num2str(IFconvtol)]);

%--------------------------------------------------------------------------
% create filenames for input files from type '*.mat'
filename1 = fullfile(pwd, 'comp_geo', 'xfeminputdata_comp_geo.mat');
filename2 = fullfile(pwd, 'preprocess', 'xfeminputdata_preprocess.mat');
filename3 = fullfile(pwd, 'XFEM', 'xfeminputdata_xfem.mat');

% save input parameters to mat-files into specific directories
save(filename1, 'IFmeshstructure', 'IFshapegeometryID', 'IFlength', ...
    'IFheight', 'IFnldivx', 'IFnldivy', 'IFdatasetp');  % for 'comp_geo'
save(filename2, 'IFDirichletBCs', 'IFNeumannBCs', 'IFMatSet'); % for 'preprocess'
save(filename3, 'IFsliding_switch','IFmethod','IFpenalty','IFnitsche',...
    'IFSolverType','IFmaxiter','IFconvtol');                    % for 'XFEM'

% clear workspace
clear IFmeshstructure IFshapegeometryID IFlength IFheight IFnldivx ...
    IFnldivy IFdatasetp IFDirichletBCs IFNeumannBCs IFMatSet ...
    IFsliding_switch IFmethod IFpenalty IFnitsche IFSolverType ...
    IFmaxiter IFconvtol;
clear filename1 filename2 filename3;