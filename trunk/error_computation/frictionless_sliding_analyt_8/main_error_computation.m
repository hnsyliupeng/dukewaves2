% main_error_computation.m
%
% This scrips manages the computation of different error measures.
%

% Author: Matthias Mayr (05/2010)

% print a brief message into console
disp('error computation ...')

% call subroutine to compute error in displacement, measured in L2-norm
error_norm_enr;

% call subroutine to compute error in energy
% energy_norm_enr;

% call subroutine to compute error in tractions at the interface, measured
% in L2-norm
error_norm_traction;




