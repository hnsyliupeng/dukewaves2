% applybcs.m
%
% This function manages the applying of DBCs and NBCs. These are applied
% separately. The choice, which set of BCs has to be applied,
% is made in the input file by setting 'IFDirichletBCs' and IFNeumannBCs'.
% Routine applies boundary conditions based on geometries.
%
% The different sets of BCs are stored in separate skripts (m-files), which
% are called by this method, depending on unser's choice. To keep the file
% structur clear, the BC-scripts are stored in the subfolder
% 'Boundary_Conditions'. So, before running these scripts, the current
% directory is changed by a 'cd'-command.
%
% Input Parameters:
%   x               x-coordinates of all nodes
%   y               y-coordinates of all nodes
%   numnod          number of nodes in discretization
%   beam_l
%   beam_h
%   f
%
% Returned Variables:
%   force           matrix with external nodal forces in x- and y-direction
%                   for each node   
%   dispbc          knows, in which DOFs there is a DBC    
%   ubar            stores the values of the DBCs
%   dispbc2         knows, in which DOFs there is a DBC (second set of DBCs)
%   ubar2           stores the values of the DBCs (second set of DBCs)
%   dispbc3         knows, in which DOFs there is a DBC (third set of DBCs)
%   ubar3           stores the values of the DBCs (third set of DBCs)
%   num_enr_surf
%   enr_surfs
%   bc_enr
%

% Author: Matthias Mayr (04/2010)

function [force,dispbc,ubar,dispbc2,ubar2,dispbc3,ubar3,num_enr_surf,enr_surfs,bc_enr,nodeNBC,FORCE]...
    = applybcs(x,y,numnod,beam_l,beam_h,f)

% load parameters from input file 'xfeminputdata_preprocess.mat'
load xfeminputdata_preprocess.mat

bc_enr = 0;

% Initialize
num_enr_surf = 0;
enr_surfs = struct('nodes',[],'xsi',[],'coords',[],'grain',0);
num_edges = 0;

ubar = zeros(2,numnod);
dispbc = zeros(2,numnod);
force = zeros(2,numnod);

% generate default variables for a second and third set of boundary
% conditions
dispbc2 = zeros(2,numnod);
ubar2 = zeros(2,numnod); 
dispbc3 = zeros(2,numnod);
ubar3 = zeros(2,numnod); 

% load DBCs
switch IFDirichletBCs
  case 0
    filename_DBC = fullfile(pwd,'Boundary_Conditions','bc_conv4_DBC');
  case 1
    filename_DBC = fullfile(pwd,'Boundary_Conditions','bc_conv1_DBC');
  case 2
    filename_DBC = fullfile(pwd,'Boundary_Conditions','bc_conv7_DBC');
  case 3
    filename_DBC = fullfile(pwd,'Boundary_Conditions','bc_conv9_DBC');
  case 4
    filename_DBC = fullfile(pwd,'Boundary_Conditions','bc_conv10_DBC');
  case 5
    filename_DBC = fullfile(pwd,'Boundary_Conditions','bc_conv11_DBC');
  case 6
    filename_DBC = fullfile(pwd,'Boundary_Conditions','multi1_DBC');
  case 7
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding1_24_6_DBC');
  case 8
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding1_72_18_DBC');
  case 9
    filename_DBC = fullfile(pwd,'Boundary_Conditions','quarterring_gmsh_DBC');
  case 10
    filename_DBC = fullfile(pwd,'Boundary_Conditions','patchtest_8_2_DBC');
  case 11
    filename_DBC = fullfile(pwd,'Boundary_Conditions','patchtest_24_6_DBC');
  case 12
    filename_DBC = fullfile(pwd,'Boundary_Conditions','patchtest_72_18_DBC');
  case 13
    filename_DBC = fullfile(pwd,'Boundary_Conditions','patchtest_4_1_DBC');
  case 14
    filename_DBC = fullfile(pwd,'Boundary_Conditions','patchtest_392_98_DBC');
  case 15
    filename_DBC = fullfile(pwd,'Boundary_Conditions','bc_conv7_frictionless_sliding_DBC');
  case 16
    filename_DBC = fullfile(pwd,'Boundary_Conditions','Hertzian_Contact_1_DBC');
  case 17
    filename_DBC = fullfile(pwd,'Boundary_Conditions','Hertzian_Contact_2_DBC');
  case 18
    filename_DBC = fullfile(pwd,'Boundary_Conditions','square_DBC');
  case 19
    filename_DBC = fullfile(pwd,'Boundary_Conditions','patchtest_gmsh_1604_DBC');
  case 20
    filename_DBC = fullfile(pwd,'Boundary_Conditions','patchtest_gmsh_130_DBC');
  case 21
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_1_130_DBC');
  case 22
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_1_1604_DBC');
  case 23
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_1_46_DBC');
  case 24
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_1_534_DBC');
  case 25
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_1_3622_DBC');
  case 26
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_1_6660_DBC');
  case 27
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_1_72_18_DBC');
  case 28
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_2_552_DBC');
  case 29
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_3_1053_DBC');
  case 30
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_4_328_DBC');
  case 31
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_4_1840_DBC');
  case 32
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_4_29229_DBC');
  case 33
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_5_144_DBC');
  case 34
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_6_81_6_DBC');
  case 35
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_7_DBC');
  case 36
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_4_7337_DBC');
  case 37
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_8_81_12_DBC');
  case 38
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_2_300_DBC');
  case 39
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_2_9188_DBC');
  case 40
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_2_72_18_DBC');
  case 41
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_2_122_DBC');
  case 42
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_2_246_DBC');
  case 43
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_2_814_DBC');
  case 44
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_2_1854_DBC');
  case 45
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_2_4758_DBC');
  case 46
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_2_14398_DBC');
  case 47
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_2_6220_DBC');
  case 48
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_2_31424_DBC');
  case 49
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_6_41_3_DBC');
  case 50
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_6_161_12_DBC');
  case 51
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_6_321_24_DBC');
  case 52
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_6_641_48_DBC');
  case 53
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_6_1281_96_DBC');
  case 54
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_8_534_DBC');
  case 55
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_8_3622_DBC');
  case 56
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_8_41_6_DBC');
  case 57
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_8_81_20_DBC');
  case 58
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_8_121_30_DBC');
  case 59
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_8_161_40_DBC');
  case 60
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_8_321_80_DBC');
  case 61
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_8_41_10_DBC');
  case 62
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_8_21_6_DBC');
  case 63
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_8_1604_DBC');
  case 64
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_8_6660_DBC');
  case 65
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_8_23174_DBC');
  case 66
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding1_6660_DBC');
  case 67
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding1_2020_DBC');
  case 68
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_8_420_DBC');
  case 69
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_8_3706_DBC');
  case 70
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_8_14766_DBC');
  case 71
    filename_DBC = fullfile(pwd,'Boundary_Conditions','plasticity_1_7_7_DBC');
  case 72
    filename_DBC = fullfile(pwd,'Boundary_Conditions','plasticity_2_21_10_DBC');
  case 73
    filename_DBC = fullfile(pwd,'Boundary_Conditions','plasticity_2_ref_50_50_DBC');
  case 74
    filename_DBC = fullfile(pwd,'Boundary_Conditions','plasticity_2_41_20_DBC');
  case 75
    filename_DBC = fullfile(pwd,'Boundary_Conditions','plasticity_2_81_40_DBC');
  case 76
    filename_DBC = fullfile(pwd,'Boundary_Conditions','plasticity_3_125_5_DBC');
  case 77
    filename_DBC = fullfile(pwd,'Boundary_Conditions','plasticity_3_250_11_DBC');
  case 78
    filename_DBC = fullfile(pwd,'Boundary_Conditions','plasticity_3_500_21_DBC');
  case 79
    filename_DBC = fullfile(pwd,'Boundary_Conditions','plasticity_3_750_31_DBC');
  case 80
    filename_DBC = fullfile(pwd,'Boundary_Conditions','plasticity_3_1000_41_DBC');
  case 81
    filename_DBC = fullfile(pwd,'Boundary_Conditions','InputFileRoutine_DBC');
  case 82
    filename_DBC = fullfile(pwd,'Boundary_Conditions','plasticity_4_40_21_DBC');
  case 83
    filename_DBC = fullfile(pwd,'Boundary_Conditions','plasticity_4_80_41_DBC');
  case 84
    filename_DBC = fullfile(pwd,'Boundary_Conditions','Simone2006_147_49_DBC');
  case 85
    filename_DBC = fullfile(pwd,'Boundary_Conditions','Simone2006_14406_DBC');
  case 86
    filename_DBC = fullfile(pwd,'Boundary_Conditions','plasticity_5_21_27_DBC');
  case 87
    filename_DBC = fullfile(pwd,'Boundary_Conditions','plasticity_5_61_77_DBC');
  case 88
    filename_DBC = fullfile(pwd,'Boundary_Conditions','plasticity_5_ref_21_27_DBC');
  case 89
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_8_3598_DBC');
  case 90
    filename_DBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_8_104_DBC');
  otherwise
    error('MATLAB:preprocess:applybcs','Unvalid ID for Dirichlet BCs. Either change ID in input file or introduce additional case in "applybcs.m"');
end;
run(filename_DBC);

% load NBCs
switch IFNeumannBCs
  case 0
    filename_NBC = fullfile(pwd,'Boundary_Conditions','bc_conv4_NBC');
  case 1
    filename_NBC = fullfile(pwd,'Boundary_Conditions','bc_conv1_NBC');
  case 2
    filename_NBC = fullfile(pwd,'Boundary_Conditions','bc_conv7_NBC');
  case 3
    filename_NBC = fullfile(pwd,'Boundary_Conditions','bc_conv9_NBC');
  case 4
    filename_NBC = fullfile(pwd,'Boundary_Conditions','bc_conv10_NBC');
  case 5
    filename_NBC = fullfile(pwd,'Boundary_Conditions','bc_conv11_NBC');
  case 6
    filename_NBC = fullfile(pwd,'Boundary_Conditions','multi1_NBC');
  case 7
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding1_24_6_NBC');
  case 8
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding1_72_18_NBC');
  case 9
    filename_NBC = fullfile(pwd,'Boundary_Conditions','quarterring_gmsh_NBC');
  case 10
    filename_NBC = fullfile(pwd,'Boundary_Conditions','patchtest_8_2_NBC');
  case 11
    filename_NBC = fullfile(pwd,'Boundary_Conditions','patchtest_24_6_NBC');
  case 12
    filename_NBC = fullfile(pwd,'Boundary_Conditions','patchtest_72_18_NBC');
  case 13
    filename_NBC = fullfile(pwd,'Boundary_Conditions','patchtest_4_1_NBC');
  case 14
    filename_NBC = fullfile(pwd,'Boundary_Conditions','patchtest_392_98_NBC');
  case 15
    filename_NBC = fullfile(pwd,'Boundary_Conditions','hatshaped_x_forces_72_18_NBC');
  case 16
    filename_NBC = fullfile(pwd,'Boundary_Conditions','parabolic_x_forces_72_18_NBC');
  case 17
    filename_NBC = fullfile(pwd,'Boundary_Conditions','Hertzian_Contact_1_NBC');
  case 18
    filename_NBC = fullfile(pwd,'Boundary_Conditions','patchtest_enr_NBC_8_2_NBC');   
  case 19
    filename_NBC = fullfile(pwd,'Boundary_Conditions','Hertzian_Contact_2_NBC');
  case 20
    filename_NBC = fullfile(pwd,'Boundary_Conditions','square_NBC');
  case 21
    filename_NBC = fullfile(pwd,'Boundary_Conditions','patchtest_gmsh_1604_NBC');
  case 22
    filename_NBC = fullfile(pwd,'Boundary_Conditions','patchtest_gmsh_130_NBC');
  case 23
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_1_130_NBC');
  case 24
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_1_1604_NBC');
  case 25
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_1_46_NBC');
  case 26
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_1_534_NBC');
  case 27
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_1_3622_NBC');
  case 28
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_1_6660_NBC');
  case 29
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_1_72_18_NBC');
  case 30
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_2_552_NBC');
  case 31
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_3_1053_NBC');
  case 32
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_4_328_NBC');
  case 33
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_4_1840_NBC');
  case 34
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_4_29229_NBC');
  case 35
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_5_144_NBC');
  case 36
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_6_81_6_NBC');
  case 37
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_7_NBC');
  case 38
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_4_7337_NBC');
  case 39
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_8_81_12_NBC');
  case 40
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_2_300_NBC');
  case 41
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_2_9188_NBC');
  case 42
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_2_72_18_NBC');
  case 43
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_2_122_NBC');
  case 44
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_2_246_NBC');
  case 45
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_2_814_NBC');
  case 46
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_2_1854_NBC');
  case 47
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_2_4758_NBC');
  case 48
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_2_14398_NBC');
  case 49
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_2_6220_NBC');
  case 50
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_2_31424_NBC');
  case 51
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_6_41_3_NBC');
  case 52
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_6_161_12_NBC');
  case 53
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_6_321_24_NBC');
  case 54
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_6_641_48_NBC');
  case 55
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_6_1281_96_NBC');
  case 56
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_8_534_NBC');
  case 57
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_8_3622_NBC');
  case 58
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_8_41_6_NBC');
  case 59
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_8_81_20_NBC');
  case 60
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_8_121_30_NBC');
  case 61
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_8_161_40_NBC');
  case 62
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_8_321_80_NBC');
  case 63
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_8_41_10_NBC');
  case 64
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_8_21_6_NBC');
  case 65
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_8_1604_NBC');
  case 66
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_8_6660_NBC');
  case 67
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_8_23174_NBC');
  case 68
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding1_6660_NBC');
  case 69
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding1_2020_NBC');
  case 70
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_8_420_NBC');
  case 71
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_8_3706_NBC');
  case 72
    filename_NBC = fullfile(pwd,'Boundary_Conditions','frictionless_sliding_analyt_8_14766_NBC');
  case 73
    filename_NBC = fullfile(pwd,'Boundary_Conditions','plasticity_1_7_7_NBC');
  case 74
    filename_NBC = fullfile(pwd,'Boundary_Conditions','no_external_loads');
  case 75
    filename_NBC = fullfile(pwd,'Boundary_Conditions','plasticity_2_21_10_NBC');
  case 76
    filename_NBC = fullfile(pwd,'Boundary_Conditions','plasticity_2_ref_50_50_NBC');
  case 77
    filename_NBC = fullfile(pwd,'Boundary_Conditions','plasticity_2_41_20_NBC');
  case 78
    filename_NBC = fullfile(pwd,'Boundary_Conditions','plasticity_2_81_40_NBC');
  case 79
    filename_NBC = fullfile(pwd,'Boundary_Conditions','plasticity_3_125_5_NBC');
  case 80
    filename_NBC = fullfile(pwd,'Boundary_Conditions','plasticity_3_250_11_NBC');
  case 81
    filename_NBC = fullfile(pwd,'Boundary_Conditions','plasticity_3_500_21_NBC');
  case 82
    filename_NBC = fullfile(pwd,'Boundary_Conditions','plasticity_3_750_31_NBC');
  case 83
    filename_NBC = fullfile(pwd,'Boundary_Conditions','plasticity_3_1000_41_NBC');
  case 84
    filename_NBC = fullfile(pwd,'Boundary_Conditions','InputFileRoutine_NBC');
  case 85
    filename_NBC = fullfile(pwd,'Boundary_Conditions','plasticity_4_40_21_NBC');
  case 86
    filename_NBC = fullfile(pwd,'Boundary_Conditions','Simone2006_147_49_NBC');
  case 87
    filename_NBC = fullfile(pwd,'Boundary_Conditions','Simone2006_14406_NBC');
  case 88
    filename_NBC = fullfile(pwd,'Boundary_Conditions','plasticity_5_ref_21_27_NBC');
  case 89
    filename_NBC = fullfile(pwd,'Boundary_Conditions','plasticity_5_61_77_NBC');
  case 90
    filename_NBC = fullfile(pwd,'Boundary_Conditions','plasticity_4_80_41_NBC');
  otherwise
    error('MATLAB:preprocess:UnvalidID', ...
      'Unvalid ID for Neumann BCs. Either change ID in input file or introduce additional case in "applybcs.m"');
end;
run(filename_NBC);

if exist('nodeNBC','var') == 0
  nodeNBC = zeros(1,3);
end;

        % set default values, which are never used, because NBCs are given nodally.
if exist('FORCE','var') == 0
    FORCE = 0;
end;
if exist('nodeNBC','var') == 0
    nodeNBC = 0;
end;
  
% clear some variables from workspace
clear IFDirichletBCs IFNeumannBCs;

