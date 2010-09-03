% getlsparameter.m
%
% CALL: getlsparameter(stol,residual_old,deltanewton, ...
%     bigk_el,seg_cut_info,INTERFACE_MAP,x,y, ...
%     big_force_loadstep,node,totaldis,IFpenalty_normal, ...
%     IFpenalty_tangential,IFmethod,IFsliding_switch,IFintegral, ...
%     IFyieldstress,id_eqns,id_dof,deltaload,DBCmatrix,IFnitsche_normal,...
%       IFnitsche_tangential,dis,old_ndisp,cutlist,maxngrains, ...
%       GRAININFO_ARR,nodegrainmap,youngs,poissons,dis_conv, ...
%       old_ndisp_conv,IFsymmetrized)
%
% Compute the line search parameter according to the algorithm proposed in
% 'Matthies, H and Strang, G: The Solution Of Nonlinear Finite Element 
% Equations (1979)'. Notations and variable names are adopted from that
% paper.
%
% Input parameters:
%   stol                  tolerance for line search
%   residual_old          residual before updating the displacements
%   deltanewton           displacement increment of current Newton step
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
%   DBCmatrix             information about Dirichlet boundary conditions
% 
% Returned variables:
%   ls_par          line search parameter
%

% Author: Matthias Mayr (08/2010)

function [ls_par] = getlsparameter(stol,residual_old,deltanewton, ...
    bigk_el,seg_cut_info,INTERFACE_MAP,x,y, ...
    big_force_loadstep,node,totaldis,IFpenalty_normal, ...
    IFpenalty_tangential,IFmethod,IFsliding_switch,IFintegral, ...
    IFyieldstress,id_eqns,id_dof,deltaload,DBCmatrix,IFnitsche_normal,...
      IFnitsche_tangential,dis,old_ndisp,cutlist,maxngrains, ...
      GRAININFO_ARR,nodegrainmap,youngs,poissons,dis_conv, ...
      old_ndisp_conv,IFsymmetrized)

% update residual
[residual ~] = ...
      buildresidual(bigk_el,seg_cut_info,INTERFACE_MAP,x,y, ...
      big_force_loadstep,node,totaldis+deltanewton,IFpenalty_normal, ...
      IFpenalty_tangential,IFmethod,IFsliding_switch,IFintegral, ...
      IFyieldstress,id_eqns,id_dof,deltaload,IFnitsche_normal,...
      IFnitsche_tangential,dis,old_ndisp,cutlist,maxngrains, ...
      GRAININFO_ARR,nodegrainmap,youngs,poissons,dis_conv, ...
      old_ndisp_conv,IFsymmetrized);
 
% impose displacement boundary conditions
% loop over all DOFs
numdof = length(totaldis);
for m=1:numdof
  % check, if current DOF has a DBC
  if DBCmatrix(1,m) == 1
    residual = residual - bigk_el(:,m) * (totaldis(m) - DBCmatrix(2,m));
    residual(m) = totaldis(m) - DBCmatrix(2,m);%totaldisDBC(m); 
  end;
end;
  
a = 0;
b = 1;
s = 1.0;

% check, if line search is necessary
Gb = deltanewton' * residual;
Ga = deltanewton' * residual_old;

if Gb * Ga < 0
%   side = 0;
  for n=1:10
    s = (Gb*a - Ga*b) / (Gb - Ga);
    if abs(b-a) < 0.0001 || abs(Gb- Ga) < 0.0001
      break;
    end;
    
    [resid ~] = ...
      buildresidual(bigk_el,seg_cut_info,INTERFACE_MAP,x,y, ...
      big_force_loadstep,node,totaldis + s * deltanewton,IFpenalty_normal, ...
      IFpenalty_tangential,IFmethod,IFsliding_switch,IFintegral, ...
      IFyieldstress,id_eqns,id_dof,deltaload,IFnitsche_normal,...
      IFnitsche_tangential,dis,old_ndisp,cutlist,maxngrains, ...
      GRAININFO_ARR,nodegrainmap,youngs,poissons,dis_conv, ...
      old_ndisp_conv,IFsymmetrized);
    
    % impose displacement boundary conditions
    % loop over all DOFs
    numdof = length(totaldis);
    for m=1:numdof
      % check, if current DOF has a DBC
      if DBCmatrix(1,m) == 1
        resid = resid - bigk_el(:,m) * (totaldis(m) + s*deltanewton(m)- DBCmatrix(2,m));
        resid(m) = totaldis(m) + s*deltanewton(m)- DBCmatrix(2,m);%totaldisDBC(m); 
      end;
    end;
    
    Gs = deltanewton' * resid;
    
    if Gs * Ga > 0
      a = s;
      Gb = Gs;
    elseif Gb * Gs > 0
      b=s;
      Gb = Gs;
    else
      break;
    end;
  end;
end

ls_par = s;

% check the admissibility of 'ls_par'
if ls_par <= 0 || ls_par > 1.0
  error('MATLAB:XFEM:Unvalidlinesearch', ...
    'Line search parameter out of admissible intervall [0;1]');
end;

%{
if abs(G) > stol * abs(G0)% || G * G0 < 0
  % initialize
  linmax = 10;
  smax = 16;
  sb = 0;
  sa = 1;
  Gb = G0;
  Ga = G;
  iteration = 0;

  % while-loop to find bracket on zero
  while ((sign(Ga) * sign(Gb)) > 0) && (sa < smax)
    sb = sa;
    sa = 2 * sa;
    Gb = Ga;
   [resid ~] = ...
      buildresidual(bigk_el,seg_cut_info,INTERFACE_MAP,x,y, ...
      big_force_loadstep,node,totaldis + sa * deltanewton,IFpenalty_normal, ...
      IFpenalty_tangential,IFmethod,IFsliding_switch,IFintegral, ...
      IFyieldstress,id_eqns,id_dof,deltaload);
    
    % impose displacement boundary conditions
    % loop over all DOFs
    numdof = length(totaldis);
    for m=1:numdof
      % check, if current DOF has a DBC
      if DBCmatrix(1,m) == 1
        resid = resid - bigk_el(:,m) * (totaldis(m) - DBCmatrix(2,m));
        resid(m) = totaldis(m) - DBCmatrix(2,m);%totaldisDBC(m); 
      end;
    end;
    
    Ga = deltanewton' * resid;
  end;

  step = sa;
  G = Ga;

  % Illinois-algorithm to find zero
  while ((sign(Ga) * sign(Gb)) > 0) && ...
    ((abs(G) > stol * abs(G0)) || (abs(sb - sa) > stol * 0.5 * (sb + sa))) && ...
    iteration < linmax

    step = sa - Ga * (sa - sb) / (Ga - Gb);

    [resid ~] = ...
        buildresidual(bigk_el,seg_cut_info,INTERFACE_MAP,x,y, ...
        big_force_loadstep,node,totaldis + step * deltanewton,IFpenalty_normal, ...
        IFpenalty_tangential,IFmethod,IFsliding_switch,IFintegral, ...
        IFyieldstress,id_eqns,id_dof,deltaload);
      
    % impose displacement boundary conditions
    % loop over all DOFs
    numdof = length(totaldis);
    for m=1:numdof
      % check, if current DOF has a DBC
      if DBCmatrix(1,m) == 1
        resid = resid - bigk_el(:,m) * (totaldis(m) - DBCmatrix(2,m));
        resid(m) = totaldis(m) - DBCmatrix(2,m);%totaldisDBC(m); 
      end;
    end;

    G = deltanewton' * resid;

    if sign(G) * sign(Ga) > 0
      Gb = 0.5 * Gb;
    else
      sb = sa;
      Gb = Ga;
    end;

    sa = step;
    Ga = G;
    
    iteration = iteration + 1;
  end;
  
  ls_par = sa;
else
  ls_par = 1;
end;
%}
end

