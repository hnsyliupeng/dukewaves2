% buildresidual.m
%
% CALL: buildresidual(bigk_el,seg_cut_info,INTERFACE_MAP,x,y, ...
%         big_force_loadstep,node,totaldis,IFpenalty_normal, ...
%         IFpenalty_tangential,IFmethod,IFsliding_switch,IFintegral, ...
%         IFyieldstress,id_eqns,id_dof,deltaload)
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
    IFnitsche_tangential,dis,orig_ndisp,cutlist,maxngrains, ...
    GRAININFO_ARR,nodegrainmap)

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
            case 0  % fully tied case
              penalty_normal = IFnitsche_normal;
              penalty_tangent = IFnitsche_tangential;
            case 1  % frictionless sliding
              penalty_normal = IFnitsche_normal;
              penalty_tangent = 0;
            case 2  % perfect plasticity
              % provide both penalty parameters: the normal one will be
              % used indepentent of the slidestate, the tangential one
              % is needed for the return mapping algorithm.
              penalty_normal = IFnitsche_normal;
              penalty_tangent = IFnitsche_tangential;
            otherwise
              error('MATLAB:XFEM','Unvalid Sliding-ID');
          end;
          
          % get penalty-residual contributions vectors, but choose between the
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

          % get Nitsche contributions
          res_nitsche = get_ele_residual_nitsche(xcoords,ycoords, ...
            seg_cut_info(i,e),INTERFACE_MAP(i).endpoints,node,x,y, ...
            dis,orig_ndisp,id_dof,cutlist,maxngrains,totaldis', ...
            id_eqns(elenodes,:),GRAININFO_ARR,nodegrainmap);          
          
          
          % compute ele_residual_constraint
          ele_residual_penalty = res_penalty_normal ...  % normal traction
                               + res_penalty_tangent; ... % tangential traction
                               
          % Since the penalty contributions act only int the enriched
          % degrees of freedom, these have to be filled with zeros
          ele_residual_constraint = [zeros(6,1); ele_residual_penalty] ... 
                                  + res_nitsche; 
                                
          % add DOFs for base DOFs to 'id'
          id_array = [zeros(1,6) id];
          id_array(1) = id_eqns(elenodes(1),1);
          id_array(2) = id_eqns(elenodes(1),2);
          id_array(3) = id_eqns(elenodes(2),1);
          id_array(4) = id_eqns(elenodes(2),2);
          id_array(5) = id_eqns(elenodes(3),1);
          id_array(6) = id_eqns(elenodes(3),2);
          
          % assemble into global residual
          for a = 1:length(id_array)
            if id_array(a) ~= 0
              residual(id_array(a)) = residual(id_array(a)) ...
                                    + ele_residual_constraint(a);
            end;
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

