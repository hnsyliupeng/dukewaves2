% returnmappingnitsche.m
%
% CALL: returnmappingnitsche(alpha_t,tgapplconv, ...
%         ttracconv,yieldstress,tgapcurrent)
%
% Input parameters:
%   alpha_t         tangential stabilizatiion parameter
%   tgapplconv      converged tangential gap of previous load step
%   ttracconv       converged tangential traction of previous load step
%   yieldstress     yield stress
%   tgapcurrent     increment in tangential gap in comparison to the last
%                   converged load step (equals the "given strain
%                   increment" in 'Simo1998'
%
% Returned variables
%   ttrac           tangential traction
%   tgappl          plastic part of tangential gap
%   f_trial         flow rule evaluated at trial state (assumed elasticity)
%

% Author: Matthias Mayr (07/2010)

function [ttrac tgappl f_trial] = returnmappingnitsche(alpha_t, ...
  tgapplconv,ttracconv,yieldstress,tgapcurrent, ...
  delta_stress_avg_tang_scalar,delta_stress_avg,normal,tangent,gap_inc, ...
  n_matrix,stress_avg,gap)
% The computation is done in a reference frame attaced to the interface, so
% it is sufficient to compute the scalar value of the tangential traction.
% The vector in the 2D-plane can be obtained by projecting the scalar value
% into the tangential direction.
%% compute the trial state
% trial traction
% ttrac_trial = ttracconv + alpha_t * tgapcurrent - delta_stress_avg_tang_scalar;
% ttrac_trial = ttracconv + tangent' * (alpha_t * gap_inc - n_matrix * delta_stress_avg);
ttrac_trial = alpha_t * (gap' * tangent - tgapplconv) - (n_matrix * stress_avg)' * tangent;

% trial plastic gap
tgappl_trial = tgapplconv;                        

% flow rule of trial state
f_trial = abs(ttrac_trial) - yieldstress;       
% ----------------------------------------------------------------------- %
%% Check the admissibility of the flowrule 'f_trial'
if f_trial <= 0
  % ==> elastic deformation
  % tangential traction (equal to trial state)
  ttrac = ttrac_trial;                            
  
  % tangential plastic gap (equal to trial state)  
  tgappl = tgappl_trial; 
else
  % ==> plastic deformation
  % algorithmic consistency parameter
  Deltagamma = f_trial / alpha_t; % is > 0 since f_trial > 0 and alpha_t > 0

  % tangential traction
  ttrac = ttrac_trial - alpha_t * Deltagamma * sign(ttrac_trial); 
  
  % tangential plastic gap
  tgappl = tgapplconv + Deltagamma * sign(ttrac_trial); 
  
  % flow rule (should be zero, now)
  f = abs(ttrac) - yieldstress;
end; 
% ----------------------------------------------------------------------- %
end

