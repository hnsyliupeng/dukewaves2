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

function [force,dispbc,ubar,num_enr_surf,enr_surfs,bc_enr]...
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
        filename_DBC = fullfile(pwd,'Boundary_Conditions','rect_fixabug_DBC');
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
        filename_NBC = fullfile(pwd,'Boundary_Conditions','rect_fixabug_NBC');
    otherwise
        error('MATLAB:preprocess:applybcs','Unvalid ID for Neumann BCs. Either change ID in input file or introduce additional case in "applybcs.m"');
end;
run(filename_NBC);

        
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


%unstruct2
%pressure3
%run(f)

% clear some variables from workspace
clear IFDirichletBCs IFNeumannBCs;

