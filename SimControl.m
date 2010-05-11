% SimControl.m
%
% This script controls the application flow. It manages the calls of
% different subscripts and the file transfer between the subdirectories.
%
% You can determine, which parts of the code will be executed. You
% have to choose a input file, too.
%

% Author: Matthias Mayr (04/2010)

% clear console, screen and workspace
clc; close all; clear all;

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% determine, which parts get executed. 0 = false, 1 = true
EXEC_comp_geo = 1;          % generate a new background mesh
EXEC_preprocess = 1;        % apply new BCs
EXEC_XFEM = 1;              % solve

% set the filename of the input file without file extension '.m', that has 
% to be used.
%
% filename_input_file = 'InputFileRoutine';   % NO FILE EXTENSION '.m'
filename_input_file = 'inp_square';   % NO FILE EXTENSION '.m'
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% The flow control starts here.
%--------------------------------------------------------------------------

% execute input parameter script
run(filename_input_file);

% mesh generation
if EXEC_comp_geo == 1
    % run mesh generation
    filename_comp_geo = fullfile(pwd,'comp_geo','main_comp_geo');
    run(filename_comp_geo);
    
    % copy the mesh-file into the preprocess-directory
    source = fullfile(pwd,'comp_geo','my_new_mesh.mat');
    destination = fullfile(pwd,'preprocess','my_new_mesh.mat');
    [status,message,messageid] = copyfile(source, destination);
    if messageid ~= ''
        error('MATLAB:SimControl','Copying the mesh-file "my_new_mesh.mat" failed.');
    end;
end;

% preprocessing
if EXEC_preprocess == 1
    % run preprocessing (apply BCs)
    filename_preprocess = fullfile('preprocess','main_preprocess');
    run(filename_preprocess);

    % copy the mesh-file with BCs into the XFEM-directory
    source = fullfile(pwd,'preprocess','my_new_mesh_with_BCs.mat');
    destination = fullfile(pwd,'XFEM','my_new_mesh_with_BCs.mat');
    [status,message,messageid] = copyfile(source, destination);
    if messageid ~= ''
        error('MATLAB:SimControl','Copying the mesh-file with BCs "my_new_mesh_with_BCs.mat" failed.');
    end;
end;

% XFEM-Computation and solution process
if EXEC_XFEM == 1
    % run solving and postprocessing
    filename_XFEM = fullfile(pwd,'XFEM','main_xfem');
    run(filename_XFEM);
end;

%--------------------------------------------------------------------------
% End of flow control
%--------------------------------------------------------------------------

% clear some temporary variables
clear EXEC_comp_geo EXEC_preprocess EXEC_XFEM messageid ...
    filename_comp_geo filename_input_file filename_preprocess ...
    filename_XFEM source destination status;

disp('SimControl finished succesfully.');

