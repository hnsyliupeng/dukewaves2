% Input File 'inp_frictionless_sliding1.m'
%
% Here, you can define all parameters to configure the simulation.
%
%**************************************************************************
% GIVE A SHORT DESCRIPTION OF THE EXAMPLE
%**************************************************************************
% frictionless sliding: rectangular domain with interface with frictionless
% sliding. Left part of the domain fixed, right part of the domain is
% pulled in x-direction. Left part with poisson = 0.0, right part with
% poisson = 0.3. Transversal deformation of right part causes jump in
% displacement field at the interface. Length x Height = 16 x 4. Elements
% 24 x 6.
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
%--------------------------------------------------------------------------
% PARAMETERS FOR INTERFACES
% Set some parameters to specify the interfaces (boundaries of the grains)
%
% Choose one of the datasets for p in 'comp_geo/vdata_multi.m'
%
IFdatasetp = 4;
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
IFDirichletBCs = 7;
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
IFNeumannBCs = 7;
%--------------------------------------------------------------------------
% MATERIAL PROPERTIES
% Set an ID 'IFMatSet' to chose a set of material properties from material 
% database 'preprocess\MaterialProperties.m'
% ID    Description
% 0     all grains with same properties (nue = 0.3, E = 1000.0)
% 1     Two grains (nue1 = 0.0, nue2 = 0.3, E1 = E2 = 1000.0)
IFMatSet = 1;
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

%--------------------------------------------------------------------------
% create filenames for input files from type '*.mat'
filename1 = fullfile(pwd, 'comp_geo', 'xfeminputdata_comp_geo.mat');
filename2 = fullfile(pwd, 'preprocess', 'xfeminputdata_preprocess.mat');
filename3 = fullfile(pwd, 'XFEM', 'xfeminputdata_xfem.mat');

% save input parameters to mat-files into specific directories
save(filename1, 'IFmeshstructure', 'IFshapegeometryID', 'IFlength', ...
    'IFheight', 'IFnldivx', 'IFnldivy', 'IFdatasetp');  % for 'comp_geo'
save(filename2, 'IFDirichletBCs', 'IFNeumannBCs', 'IFMatSet'); % for 'preprocess'
save(filename3, 'IFsliding_switch');                    % for 'XFEM'

% clear workspace
clear IFmeshstructure IFshapegeometryID IFlength IFheight IFnldivx ...
    IFnldivy IFdatasetp IFDirichletBCs IFNeumannBCs IFMatSet IFsliding_switch;
clear filename1 filename2 filename3;