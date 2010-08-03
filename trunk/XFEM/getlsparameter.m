% getlsparameter.m
%
% CALL: getlsparameter()
%
% Compute the line search parameter according to the algorithm proposed in
% 'Matthies, H and Strang, G: The Solution Of Nonlinear Finite Element 
% Equations (1979)'. Notations and variable names are adopted from that
% paper.
%
% Input parameters:
%
%
% Returned variables:
%   ls_par          line search parameter
%

% Author: Matthias Mayr (08/2010)

function [ls_par] = getlsparameter(stol,Gb,Ga,d)

% Initilize
linmax = 10;
smax = 16;
sb = 0;
sa = 1;

% while-loop
while ((sign(Ga) * sign(Gb)) > 0) && (sa < smax)
  sb = sa;
  sa = 2 * sa;
  Gb = Ga;
  resid = ;
  Ga = d * resid;
end;

step = Ga;
G = Ga;


while ((sign(Ga) * sign(Gb)) > 0) && ...
    ((abs(G) > stol * abs(G0)) || (abs(sb - sa) > stol * 0.5 * (sb + sa)))

step = sa - Ga * (sa - sb) / (Ga - Gb)
resid = ;
G = d * resid;

if sign(G) * sign(Ga) > 0
  Gb = 0.5 * Gb;
else
  sb = sa;
  Gb = Ga;
end;

sa = step;
Ga = G;
end;

end

