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
%   x
%   y
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
%   num_enr_surf
%   enr_surfs
%   bc_enr
%

% Author: Matthias Mayr (04/2010)

function [force,dispbc,ubar,num_enr_surf,enr_surfs,bc_enr,nodeNBC,FORCE]...
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
  otherwise
    error('MATLAB:preprocess:applybcs','Unvalid ID for Neumann BCs. Either change ID in input file or introduce additional case in "applybcs.m"');
end;
run(filename_NBC);

if exist('nodeNBC','var') == 0
  nodeNBC = zeros(1,3);
end;

        
%  for i=1:numnod
%    if (x(i) == 0)
%         dispbc(1,i) = 1.0;
%         dispbc(2,i) = 1.0;
%    end
%    
%    if (x(i) == beam_l)
%         dispbc(1,i) = 1.0;
%         ubar(1,i) = 0.9;
%         if y(i) == 0
%             dispbc(2,i) = 1.0;
%             ubar(2,i) = 0.0;
%         end
%    end
% 
%  end


% set default values, which are never used, because NBCs are given nodally.
if exist('FORCE','var') == 0
    FORCE = 0;
end;
if exist('nodeNBC','var') == 0
    nodeNBC = 0;
end;

%unstruct2
%pressure3
%run(f)

% clear some variables from workspace
clear IFDirichletBCs IFNeumannBCs;

