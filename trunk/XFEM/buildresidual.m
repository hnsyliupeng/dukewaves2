% buildresidual.m
%
% CALL: buildresidual(bigk_el,seg_cut_info,INTERFACE_MAP,x,y, ...
%         big_force_loadstep,node,totaldis,IFpenalty_normal, ...
%         IFpenalty_tangential,IFmethod,IFsliding_switch,IFintegral, ...
%         IFyieldstress,id_eqns,id_dof,deltaload,IFnitsche_normal, ...
%         IFnitsche_tangential,dis,old_ndisp,cutlist,maxngrains, ...
%         GRAININFO_ARR,nodegrainmap,youngs,poissons)
%
% This subroutine assembles the global residual. All variable names are 
% adopted from 'main_xfem.m' to preserve consistency.
%
% Input parameters:
%   bigk_el               elastic contribution to global tangent matrix
%   seg_cut_info          information about cut elements
%   INTERFACE_MAP         information about interfaces
%   x                     x-coordinates of all nodes
%   y                     y-coordinates of all nodes
%   big_force_loadstep    external force vector for current load step
%   node                  connectivity between nodes and elements
%   totaldis              current global displacement vector
%   IFpenalty_normal      normal penalty parameter
%   IFpenalty_tangential  tangential penalty paramter
%   IFmethod              ID for method of constraint enforcement
%   IFsliding_switch      ID for sliding case
%   IFintegral            ID for number on integrals for penalty term
%   IFyieldstress         yield stress
%   id_eqns               mapping between nodes and DOFs
%   id_dof                mapping between nodes and enriching grains
%   deltaload             current displacement increment in comparison to
%                         previous loadstep
%   IFnitsche_normal      stabilization parameter in normal direction
%   IFnitsche_tangential  stabilization parameter in tangential direction
%   dis                   current re-assembled displacemement vector
%   old_ndisp             unmodified total solution vector
%   cutlist               list of cut elements
%   maxngrains            total number of all grains
%   GRAININFO_ARR         some data about all grains
%   nodegrainmap          mapping between nodes and grains
%   youngs                array with Young's moduli of all grains
%   poissons              array with Poisson's ratios of all grains
%   dis_conv              total displacement vector of previous converged
%                         load step
%   old_ndisp_conv        unmodified total solution vector of previous 
%                         converged load step
%   IFsymmetrized         unsymmetric or symmetrized version of Nitsche's
%                         method with plasticity
% 
% Returned variables:
%   residual              residual
%   seg_cut_info          updated information about cut elements
%

% Author: Matthias Mayr (08/2010)

function [residual seg_cut_info] = ...
  buildresidual(bigk_el,seg_cut_info,INTERFACE_MAP,x,y, ...
    big_force_loadstep,node,totaldis,IFpenalty_normal, ...
    IFpenalty_tangential,IFmethod,IFsliding_switch,IFintegral, ...
    IFyieldstress,id_eqns,id_dof,deltaload,IFnitsche_normal, ...
    IFnitsche_tangential,dis,old_ndisp,cutlist,maxngrains, ...
    GRAININFO_ARR,nodegrainmap,youngs,poissons,dis_conv,old_ndisp_conv, ...
    IFsymmetrized)

% contribution of elastic bulk field of the inner of the grains which
% will not be affected from plasticity at the interface
%   residual = bigk_el * totaldisDBC;
residual = bigk_el * totaldis;

% residual contributions due to constraints
for i=1:size(seg_cut_info,1)    % loop over all interfaces 'i'
  for e=1:size(seg_cut_info,2)  % loop over all cut elements 'e'
    if seg_cut_info(i,e).elemno ~= -1  % only, if 'e' is cut by 'i'
      % get some element data
      eleID = seg_cut_info(i,e).elemno; % global element ID
      elenodes = node(:,eleID); % global node IDs of this element's nodes
      xcoords = x(elenodes);    % x-coordinates of 'elenodes'
      ycoords = y(elenodes);    % y-coordinates of 'elenodes'

      switch IFmethod
        case 0  % Lagrange multipliers
        case 1  % Penalty method
          % get penalty parameters (depending on sliding case)
          switch IFsliding_switch
            case 0  % fully tied case
              penalty_normal = IFpenalty_normal;
              penalty_tangent = IFpenalty_tangential;
            case 1  % frictionless sliding
              penalty_normal = IFpenalty_normal;
              penalty_tangent = 0;
            case 2  % perfect plasticity
              % provide both penalty parameters: the normal one will be
              % used indepentent of the slidestate, the tangential one
              % is needed for the return mapping algorithm.
              penalty_normal = IFpenalty_normal;
              penalty_tangent = IFpenalty_tangential;
            otherwise
              error('MATLAB:XFEM','Unvalid Sliding-ID');
          end;

          % get residual contributions vectors, but choose between the
          % one- and two-integral formulation
          switch IFintegral
            case 1  % one integral
              [res_penalty_normal res_penalty_tangent id tgappl ttrac f_trial] = ...
                get_ele_residual_penalty_alternative(xcoords,ycoords, ...
                seg_cut_info(i,e),INTERFACE_MAP(i).endpoints, ...
                id_dof(elenodes,:),id_eqns(elenodes,:),totaldis, ...
                penalty_normal,penalty_tangent,IFyieldstress, ...
                deltaload,IFsliding_switch);
            case 2  % two integrals
              [res_penalty_normal res_penalty_tangent id] = ...
                get_ele_residual_penalty(xcoords,ycoords, ...
                seg_cut_info(i,e),INTERFACE_MAP(i).endpoints, ...
                id_dof(elenodes,:),id_eqns(elenodes,:),totaldis);
            otherwise
              error('MATLAB:XFEM:UnvalidID', ...
                'Unvalid number of integrals');
          end;

          % store plastic contribution to tangential gap
          seg_cut_info(i,e).tgappl = tgappl;

          % store scalar values of current tangential traction at each
          % gauss point
          seg_cut_info(i,e).ttrac = ttrac;

          % store evaluated flow rules of trial state at both gauss
          % points
          seg_cut_info(i,e).f_trial = f_trial;

          % compute ele_residual_constraint
          ele_residual_constraint = res_penalty_normal ...  % normal traction
                                  + res_penalty_tangent;    % tangential traction

          % assemble into global residual
          for a = 1:length(id)
            if id(a) ~= 0
              residual(id(a)) = residual(id(a)) + ele_residual_constraint(a);
            end;
          end;
        case 2  % Nitsche's mehtod
          % get stabilization parameters (depending on sliding case)
          switch IFsliding_switch
            % if IFnitsche > 0, then use the parameter given in the input
            % file, else compute a minimal stabilization parameter as
            % suggested in 'Dolbow2009'.
            case 0  % fully tied case
              if IFnitsche_normal >= 0
                penalty_normal = IFnitsche_normal;
              else
                penalty_normal = minstabiparameter(xcoords,ycoords, ...
                  seg_cut_info(i,e),youngs(seg_cut_info(i,e).grains), ...
                  poissons(seg_cut_info(i,e).grains), ...
                  INTERFACE_MAP(i).endpoints,nodegrainmap(elenodes), ...
                  IFintegral);
              end;
              if IFnitsche_tangential >= 0
                penalty_tangent = IFnitsche_tangential;
              else
                penalty_tangent = minstabiparameter(xcoords,ycoords, ...
                  seg_cut_info(i,e),youngs(seg_cut_info(i,e).grains), ...
                  poissons(seg_cut_info(i,e).grains), ...
                  INTERFACE_MAP(i).endpoints,nodegrainmap(elenodes), ...
                  IFintegral);
              end;
            case 1  % frictionless sliding
              if IFnitsche_normal >= 0
                penalty_normal = IFnitsche_normal;
              else
                penalty_normal = minstabiparameter(xcoords,ycoords, ...
                  seg_cut_info(i,e),youngs(seg_cut_info(i,e).grains), ...
                  poissons(seg_cut_info(i,e).grains), ...
                  INTERFACE_MAP(i).endpoints,nodegrainmap(elenodes), ...
                  IFintegral);
              end;
              penalty_tangent = 0;
            case 2  % perfect plasticity
              % provide both penalty parameters: the normal one will be
              % used indepentent of the slidestate, the tangential one
              % is needed for the return mapping algorithm.
              if IFnitsche_normal >= 0
                penalty_normal = IFnitsche_normal;
              else
                penalty_normal = minstabiparameter(xcoords,ycoords, ...
                  seg_cut_info(i,e),youngs(seg_cut_info(i,e).grains), ...
                  poissons(seg_cut_info(i,e).grains), ...
                  INTERFACE_MAP(i).endpoints,nodegrainmap(elenodes), ...
                  IFintegral);
              end;
              if IFnitsche_tangential >= 0
                penalty_tangent = IFnitsche_tangential;
              else
                penalty_tangent = minstabiparameter(xcoords,ycoords, ...
                  seg_cut_info(i,e),youngs(seg_cut_info(i,e).grains), ...
                  poissons(seg_cut_info(i,e).grains), ...
                  INTERFACE_MAP(i).endpoints,nodegrainmap(elenodes), ...
                  IFintegral);
              end;
            otherwise
              error('MATLAB:XFEM','Unvalid Sliding-ID');
          end;

          if IFsliding_switch == 0 || IFsliding_switch == 1
            % get residual contributions vectors, but choose between the
            % one- and two-integral formulation
            switch IFintegral
              case 1  % one integral
                [res_penalty_normal res_penalty_tangent id tgappl ttrac f_trial] = ...
                  get_ele_residual_penalty_alternative(xcoords,ycoords, ...
                  seg_cut_info(i,e),INTERFACE_MAP(i).endpoints, ...
                  id_dof(elenodes,:),id_eqns(elenodes,:),totaldis, ...
                  penalty_normal,penalty_tangent,IFyieldstress, ...
                  deltaload,IFsliding_switch);
              case 2  % two integrals
                [res_penalty_normal res_penalty_tangent id] = ...
                  get_ele_residual_penalty(xcoords,ycoords, ...
                  seg_cut_info(i,e),INTERFACE_MAP(i).endpoints, ...
                  id_dof(elenodes,:),id_eqns(elenodes,:),totaldis);
              otherwise
                error('MATLAB:XFEM:UnvalidID', ...
                  'Unvalid number of integrals');
            end;


            % store plastic contribution to tangential gap
            seg_cut_info(i,e).tgappl = tgappl;


            % store scalar values of current tangential traction at each
            % gauss point
            seg_cut_info(i,e).ttrac = ttrac;


            % store evaluated flow rules of trial state at both gauss
            % points
            seg_cut_info(i,e).f_trial = f_trial;

            % compute 'ele_residual_stabi' due to stabilization
            ele_residual_stabi  = res_penalty_normal ...  % normal traction
                                + res_penalty_tangent;    % tangential traction


            % assemble stabilization contribution into global residual
            for a = 1:length(id)
              if id(a) ~= 0
                residual(id(a)) = residual(id(a)) + ele_residual_stabi(a);
              end;
            end;

            % get Nitsche contributions
            [res_nitsche id_array] = get_ele_residual_nitsche(xcoords, ...
              ycoords,seg_cut_info(i,e),INTERFACE_MAP(i).endpoints,node, ...
              x,y,dis,old_ndisp,id_dof,cutlist,maxngrains,totaldis', ...
              id_eqns(elenodes,:),GRAININFO_ARR,nodegrainmap, ...
              IFsliding_switch);          

            if length(id_array) ~= 18
              error('MATLAB:XFEM','ID-array for Nitsche residual too short.');
            end;

            % assemble Nitsche contribution into global residual
            for a = 1:length(id_array)
              if id_array(a) ~= 0
                residual(id_array(a)) = residual(id_array(a)) ...
                                      + res_nitsche(a);
              end;
            end;
          else
            % get residual contributions vectors for normal direction, but 
            % choose between the one- and two-integral formulation
            switch IFintegral
              case 1  % one integral
                [res_penalty_normal res_penalty_tangent id tgappl ttrac f_trial] = ...
                  get_ele_residual_penalty_alternative(xcoords,ycoords, ...
                  seg_cut_info(i,e),INTERFACE_MAP(i).endpoints, ...
                  id_dof(elenodes,:),id_eqns(elenodes,:),totaldis, ...
                  penalty_normal,0,IFyieldstress, ...
                  deltaload,1);
              case 2  % two integrals
                [res_penalty_normal res_penalty_tangent id] = ...
                  get_ele_residual_penalty(xcoords,ycoords, ...
                  seg_cut_info(i,e),INTERFACE_MAP(i).endpoints, ...
                  id_dof(elenodes,:),id_eqns(elenodes,:),totaldis);
              otherwise
                error('MATLAB:XFEM:UnvalidID', ...
                  'Unvalid number of integrals');
            end;

            % compute 'ele_residual_stabi' due to stabilization
            ele_residual_stabi  = res_penalty_normal; ...  % normal traction

            % assemble stabilization contribution into global residual
            for a = 1:length(id)
              if id(a) ~= 0
                residual(id(a)) = residual(id(a)) + ele_residual_stabi(a);
              end;
            end;

            % get Nitsche contributions for normal direction
            [res_nitsche id_array] = get_ele_residual_nitsche(xcoords, ...
              ycoords,seg_cut_info(i,e),INTERFACE_MAP(i).endpoints,node, ...
              x,y,dis,old_ndisp,id_dof,cutlist,maxngrains,totaldis', ...
              id_eqns(elenodes,:),GRAININFO_ARR,nodegrainmap, ...
              1);          

            if length(id_array) ~= 18
              error('MATLAB:XFEM','ID-array for Nitsche residual too short.');
            end;

            % assemble Nitsche contribution into global residual
            for a = 1:length(id_array)
              if id_array(a) ~= 0
                residual(id_array(a)) = residual(id_array(a)) ...
                                      + res_nitsche(a);
              end;
            end;
            
            % get residual contribution for tangential direction
            [res_nit_tang id_tang ttrac tgappl f_trial] = ...
              get_ele_residual_nitsche_plasticity(xcoords, ...
              ycoords,seg_cut_info(i,e),INTERFACE_MAP(i).endpoints,node, ...
              x,y,dis,old_ndisp,id_dof,cutlist,maxngrains,totaldis', ...
              id_eqns(elenodes,:),GRAININFO_ARR,nodegrainmap, ...
              IFsliding_switch,penalty_tangent,IFyieldstress,deltaload, ...
              dis_conv,old_ndisp_conv,IFsymmetrized); 

            % assemble tangential Nitsche contribution into global residual
            for a = 1:18
              if id_tang(a) ~= 0
                residual(id_tang(a)) = residual(id_tang(a)) + res_nit_tang(a);
              end;
            end;
            
            % store plastic contribution to tangential gap
            seg_cut_info(i,e).tgappl = tgappl;

            % store scalar values of current tangential traction at each
            % gauss point
            seg_cut_info(i,e).ttrac = ttrac;

            % store evaluated flow rules of trial state at both gauss
            % points
            seg_cut_info(i,e).f_trial = f_trial;     
          end;     
        otherwise
          error('MATLAB:XFEM','Unvalid Method-ID');
      end;
    end;
  end;
end;

% subtract the external force vector
residual = residual - big_force_loadstep;

end

