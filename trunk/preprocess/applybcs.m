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
force = zeros(2,numnod);
dispbc = zeros(2,numnod);

% load DBCs
switch IFDirichletBCs
    case 0
        cd Boundary_Conditions
        bc_conv4_DBC
        cd ..
    case 1
    case 2
    otherwise
        error('MATLAB:preprocess:applybcs','Unvalid ID for Dirichlet BCs. Either change ID in input file or introduce additional case in "applybcs.m"');
end;

% load NBCs
switch IFNeumannBCs
    case 0
        cd Boundary_Conditions
        bc_conv4_NBC
        cd ..
    case 1
    case 2
    otherwise
        error('MATLAB:preprocess:applybcs','Unvalid ID for Neumann BCs. Either change ID in input file or introduce additional case in "applybcs.m"');
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


%unstruct2
%pressure3
%run(f)

% clear some variables from workspace
clear IFDirichletBCs IFNeumannBCs;

