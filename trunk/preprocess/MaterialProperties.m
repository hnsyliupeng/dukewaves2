% MaterialProperties.m
%
% This file is like a database for material properties. There are different
% sets of material properties. To apply one of these, you have to give its
% number in the input file parameter 'IFMatSet'.
%

% Author: Matthias Mayr (04/2010)

function[MatProp]=MaterialProperties(IFMatSet)

% Apply material parameters
switch IFMatSet
   case 0
        poissons = [0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 ...
            0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3];           % poissons rates
        youngs = [1000.0 1000.0 1000.0 1000.0 1000.0 1000.0 1000.0 1000.0...
            1000.0 1000.0 1000.0 1000.0 1000.0 1000.0 1000.0 1000.0 ...
            1000.0 1000.0 1000.0 1000.0 1000.0 1000.0 1000.0 1000.0];   % Young's moduli
    otherwise
        error('MATLAB:preprocess:MaterialProperties',...
            'Unvalid Material ID "IFMatSet". Check ID or add additional case to switch-case-structure in "MaterialProperties.m".');
end;

MatProp.poissons = poissons;
MatProp.youngs = youngs;

clear poissons youngs;

end