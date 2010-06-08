% Input File 'inp_frictionless_sliding_analyt_8_161_40.m'
%
% Here, you can define all parameters to configure the simulation.
%
%**************************************************************************
% GIVE A SHORT DESCRIPTION OF THE EXAMPLE
%**************************************************************************
% Example to compare with analytical solution to show oszillations and
% their stabilization for frictionless sliding. Rectangular domain. 
% Length x height = 16 x 4. structured mesh. 161 x 40 elements. Pure
% bending.
%**************************************************************************
%
% To set up a new example, build it in this file, so that all IDs are
% strored in this file. Save it to an own file, afterwards.
%

% Author: Matthias Mayr (05/2010)

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
IFnldivx = 161;
IFnldivy = 40;
%
% filename for boundary description file for structured meshing and NBCs 
% via integration
IFboundarydescription = 'rectangular_domain_BDes'; %NO FILE EXTENSION '.m'
%
% filename of msh-file withput file extension '.msh'
% (if reading mesh from gmsh-msh-file)
IFfilename_msh_file = 'patchtest_3622';      % NO FILE EXTENSION '.msh'
%--------------------------------------------------------------------------
% PARAMETERS FOR INTERFACES
% Set some parameters to specify the interfaces (boundaries of the grains)
%
% Choose one of the datasets for p in 'comp_geo/vdata_multi.m'
%
IFdatasetp = 25;%25;%19;
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
% 15    bc_conv7_frictionless_sliding_DBC.m
% 16    Hertzian_Contact_1_DBC.m
% 17    Hertzian_Contact_2_DBC.m
% 18    square_DBC.m
% 19    patchtest_gmsh_1604_DBC.m
% 20    patchtest_gmsh_130_DBC.m
% 21    frictionless_sliding_analyt_1_130_DBC.m
% 22    frictionless_sliding_analyt_1_1604_DBC.m
% 23    frictionless_sliding_analyt_1_46_DBC.m
% 24    frictionless_sliding_analyt_1_534_DBC.m
% 25    frictionless_sliding_analyt_1_3622_DBC.m
% 26    frictionless_sliding_analyt_1_6660_DBC.m
% 27    frictionless_sliding_analyt_1_72_18_DBC.m
% 28    frictionless_sliding_analyt_2_552_DBC.m
% 29    frictionless_sliding_analyt_3_1053_DBC.m
% 30    frictionless_sliding_analyt_4_328_DBC.m
% 31    frictionless_sliding_analyt_4_1840_DBC.m
% 32    frictionless_sliding_analyt_4_29229_DBC.m
% 33    frictionless_sliding_analyt_5_144_DBC.m
% 34    frictionless_sliding_analyt_6_80_5_DBC.m
% 35    frictionless_sliding_analyt_7_DBC.m
% 36    frictionless_sliding_analyt_4_7337_DBC.m
% 37    frictionless_sliding_analyt_8_81_13_DBC.m
% 38    frictionless_sliding_analyt_2_300_DBC.m
% 39    frictionless_sliding_analyt_2_9188_DBC.m
% 40    frictionless_sliding_analyt_2_72_18_DBC.m
% 41    frictionless_sliding_analyt_2_122_DBC.m
% 42    frictionless_sliding_analyt_2_246_DBC.m
% 43    frictionless_sliding_analyt_2_814_DBC.m
% 44    frictionless_sliding_analyt_2_1854_DBC.m
% 45    frictionless_sliding_analyt_2_4758_DBC.m
% 46    frictionless_sliding_analyt_2_14398_DBC.m
% 47    frictionless_sliding_analyt_2_6220_DBC.m
% 48    frictionless_sliding_analyt_2_31424_DBC.m
% 49    frictionless_sliding_analyt_6_41_3_DBC.m
% 50    frictionless_sliding_analyt_6_161_12_DBC.m
% 51    frictionless_sliding_analyt_6_321_24_DBC.m
% 52    frictionless_sliding_analyt_6_641_48_DBC.m
% 53    frictionless_sliding_analyt_6_1281_96_DBC.m
% 54    frictionless_sliding_analyt_8_534_DBC.m
% 55    frictionless_sliding_analyt_8_3622_DBC.m
% 56    frictionless_sliding_analyt_8_41_6_DBC.m
% 57    frictionless_sliding_analyt_8_81_20_DBC.m
% 58    frictionless_sliding_analyt_8_121_30_DBC.m
% 59    frictionless_sliding_analyt_8_161_40_DBC.m
IFDirichletBCs = 59;
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
% 17    Hertzian_Contact_1_NBC.m
% 18    patchtest_enr_NBC_8_2_NBC.m
% 19    Hertzian_Contact_2_NBC.m
% 20    square_NBC.m
% 21    patchtest_gmsh_1604_NBC.m
% 22    patchtest_gmsh_130_NBC.m
% 23    frictionless_sliding_analyt_1_130_NBC.m
% 24    frictionless_sliding_analyt_1_1604_NBC.m
% 25    frictionless_sliding_analyt_1_46_NBC.m
% 26    frictionless_sliding_analyt_1_534_NBC.m
% 27    frictionless_sliding_analyt_1_3622_NBC.m
% 28    frictionless_sliding_analyt_1_6660_NBC.m
% 29    frictionless_sliding_analyt_1_72_18_NBC.m
% 30    frictionless_sliding_analyt_2_552_NBC.m
% 31    frictionless_sliding_analyt_3_1053_NBC.m
% 32    frictionless_sliding_analyt_4_328_NBC.m% 32
% 33    frictionless_sliding_analyt_4_1840_NBC.m
% 34    frictionless_sliding_analyt_4_29229_NBC.m
% 35    frictionless_sliding_analyt_5_144_NBC.m
% 36    frictionless_sliding_analyt_6_80_5_NBC.m
% 37    frictionless_sliding_analyt_7_NBC.m
% 38    frictionless_sliding_analyt_4_7337_NBC.m
% 39    frictionless_sliding_analyt_8_81_13_NBC.m
% 40    frictionless_sliding_analyt_2_300_NBC.m
% 41    frictionless_sliding_analyt_2_9188_NBC.m
% 42    frictionless_sliding_analyt_2_72_18_NBC.m
% 43    frictionless_sliding_analyt_2_122_NBC.m
% 44    frictionless_sliding_analyt_2_246_NBC.m
% 45    frictionless_sliding_analyt_2_814_NBC.m
% 46    frictionless_sliding_analyt_2_1854_NBC.m
% 47    frictionless_sliding_analyt_2_4758_NBC.m
% 48    frictionless_sliding_analyt_2_14398_NBC.m
% 49    frictionless_sliding_analyt_2_6220_NBC.m
% 50    frictionless_sliding_analyt_2_31424_NBC.m
% 51    frictionless_sliding_analyt_6_41_3_NBC.m
% 52    frictionless_sliding_analyt_6_161_12_NBC.m
% 53    frictionless_sliding_analyt_6_321_24_NBC.m
% 54    frictionless_sliding_analyt_6_641_48_NBC.m
% 55    frictionless_sliding_analyt_6_1281_96_NBC.m
% 56    frictionless_sliding_analyt_8_534_NBC.m
% 57    frictionless_sliding_analyt_8_3622_DBC.m
% 58    frictionless_sliding_analyt_8_41_6_DBC.m
% 59    frictionless_sliding_analyt_8_81_20_DBC.m
% 60    frictionless_sliding_analyt_8_121_30_DBC.m
% 61    frictionless_sliding_analyt_8_161_40_DBC.m
IFNeumannBCs = 61;
%
% method of giving NBCs
% ID    Description
% 0     nodal forces (integration done by user, only not-enriched nodes)
% 1     tractions given as functions
IFneumann = 1;
%--------------------------------------------------------------------------
% MATERIAL PROPERTIES
% Set an ID 'IFMatSet' to chose a set of material properties from material 
% database 'preprocess\MaterialProperties.m'
% ID    Description
% 0     24 grains with same material properties (nue = 0.3, E = 1000.0)
% 1     Two grains (nue1 = 0.0, nue2 = 0.3, nue3 = 0.3, Ei = 1000.0)
% 2     24 grains with different material properties
% 3     24 grains with same material properties (nue = 0.0, E = 1000.0)
% 4     3 grains with different material properties
% 5     3 grains with same material properties (nue = 0.3, E = 200.0)
% 6     Two grains (nue1 = 0.3, nue2 = nue3 = 0.0, Ei = 1000.0)
% 7     3 grains, one of them very stiff ( E --> inf )
IFMatSet = 3;
%--------------------------------------------------------------------------
% METHOD OF ENFORCING CONSTRAINTS AT THE INTERFACE
% Set an ID to choose the method, by which the constrains shall be enforced
% at the interface
% ID    Method
% 0     Lagrange Multipliers (piecewise constant)
% 1     Penalty-Method
% 2     Nitsche's Method
IFmethod = 2;
%
% Set Penalty-Parameter
IFpenalty = 8.0e+7;
%
% Nitsche Parameter
IFnitsche = 1.0e+8;
%--------------------------------------------------------------------------
% SLIDING PARAMETERS
% Set an ID to indicate, how sliding should be treaten: 'IFsliding_switch'
% ID    Description
% 0     no slidung at all (full constraint)
% 1     frictionless sliding
% 2     perfect plasticity with shear yield stress
% 3     frictional sliding with Coulomb's friction
%
IFsliding_switch = 1; 
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
