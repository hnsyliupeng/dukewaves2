% This function either sets values for boundary conditions, or calls a
% separate script to do so.  It might be worth reworking this so that those
% scripts are called directly as functions, and have all of the necessary
% information, instead of using this intermediate routine.  It's a little
% messy to have a function call a script.

function [force,dispbc,ubar,num_enr_surf,enr_surfs,bc_enr]...
    = applybcs(x,y,numnod,beam_l,beam_h,f)

bc_enr = 0;

% Initialize
num_enr_surf = 0;
enr_surfs = struct('nodes',[],'xsi',[],'coords',[],'grain',0);
num_edges = 0;

% Routine applies boundary conditions based on geometries

ubar = zeros(2,numnod);
force = zeros(2,numnod);
dispbc = zeros(2,numnod);

 for i=1:numnod
   if (x(i) == 0)
        dispbc(1,i) = 1.0;
        dispbc(2,i) = 1.0;
   end
   
   if (x(i) == beam_l)
        dispbc(1,i) = 1.0;
        ubar(1,i) = 0.9;
        if y(i) == 0
            dispbc(2,i) = 1.0;
            ubar(2,i) = 0.0;
        end
   end

 end


%unstruct2
%pressure3
%run(f)
%bc_conv4

