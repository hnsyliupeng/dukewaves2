function [force,dispbc,ubar,num_enr_surf,enr_surfs,bc_enr]...
    = applybcs(x,y,numnod,beam_l,beam_h,f)

bc_enr = 0;

% Initialize
num_enr_surf = 0;
enr_surfs = struct('nodes',[],'xsi',[],'coords',[],'grain',0);
num_edges = 0;

% % Routine applies boundary conditions based on geometries
% 
% ubar = zeros(2,numnod);
% force = zeros(2,numnod);
% dispbc = zeros(2,numnod);
% 
%  for i=1:numnod
%    if (x(i) == 0)
% %      if (x(i) == 8)
%         dispbc(1,i) = 1.0;
%         dispbc(2,i) = 1.0;
% %      end
% %      dispbc(2,i) = 1.0;
% %       if abs(y(i) - 2) < 1e-5
% %         dispbc(1,i) = 1.0;
% %       end
% %       if abs(y(i) + 2) < 1e-5
% %         dispbc(1,i) = 1.0;
% %       end
% 
%    end
%    
%    if (x(i) == 20)
%         dispbc(1,i) = 1.0;
%         ubar(1,i) = 0.05;
%         dispbc(2,i) = 1.0;
%         ubar(2,i) = 0.0;
% %        force(2,i) = -10000
%    end
% 
% %    if (x(i) == beam_l)
% %        force(1,i) = 10000/7;
% %       if (y(i) == -2.0) || (y(i) == 2.0)
% %          force(1,i) = 5000/7;
% %       end
% %    end
%  end
%  
%bc_conv1

%unstruct2
%pressure3
%run(f)
bc_conv4

