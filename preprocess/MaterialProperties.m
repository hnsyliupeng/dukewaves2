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
    case 1
        poissons = [0.0 0.3 0.0];       % poissons ratios
        youngs = [1000.0 1000.0 1000.0];  % youngs moduli
    case 2
        poissons = [0.3 0.4 0.0 0.2 0.1 0.3 0.3 0.2 0.0 0.4 0.1 0.3 ...
            0.3 0.4 0.0 0.2 0.1 0.3 0.3 0.2 0.0 0.4 0.1 0.3];           % poissons ratios
        youngs = [1000.0 1200.0 800.0 1300.0 14000.0 1000.0 700.0 1100.0...
            1000.0 1200.0 800.0 1300.0 14000.0 1000.0 700.0 1100.0 ...
            1000.0 1200.0 800.0 1300.0 14000.0 1000.0 700.0 1100.0];   % Young's moduli
    case 3
        poissons = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ...
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0];           % poissons ratios
        youngs = [1000.0 1000.0 1000.0 1000.0 1000.0 1000.0 1000.0 1000.0...
            1000.0 1000.0 1000.0 1000.0 1000.0 1000.0 1000.0 1000.0 ...
            1000.0 1000.0 1000.0 1000.0 1000.0 1000.0 1000.0 1000.0];   % Young's moduli
    case 4
        poissons = [0.28 0.3 0.32];     % poissons ratios
        youngs = [1000.0 950.0 1050.0]; % Young's moduli
    case 5
        poissons = [0.3 0.3 0.3];       % poissons ratios
        youngs = [200.0 200.0 200.0];   % Young's moduli
    otherwise
        error('MATLAB:preprocess:MaterialProperties',...
            'Unvalid Material ID "IFMatSet". Check ID or add additional case to switch-case-structure in "MaterialProperties.m".');
end;

MatProp.poissons = poissons;
MatProp.youngs = youngs;

clear poissons youngs;

end