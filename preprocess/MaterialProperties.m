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
        poissons = [0.3 0.0 0.3];       % poissons ratios
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
    case 6
        poissons = [0.3 0.0 0.0];       % poissons ratios
        youngs = [1000.0 1000.0 1000.0];  % youngs moduli
    case 7
        poissons = [0.0 0.0 0.0];       % poissons ratios
        youngs = [1000.0 1.0e+3 1000.0];  % youngs moduli
    case 8                                                  % Paper "Chen2005"
        poissons = [0.3 0.3 0.3];       % poissons ratios
        youngs = [200e+3 200e+3 200e+3];  % youngs moduli
    case 9                                                  % Paper "Konyukhov2006"
        poissons = [0.3 0.3 0.3];       % poissons ratios
        youngs = [2.1e+4 2.1e+8 2.1e+4];  % youngs moduli
    case 10                                                 % Paper "Simone2006"
        poissons = [0.3 0.3 0.3 0.3];       % poissons ratios   
        youngs = [1.0e+4 1.0e+4 1.0e+4 1.0e+4];  % youngs moduli
    case 11
        poissons = [0.0 0.0 0.0];           % poissons ratios
        youngs = [1000.0 1000.0 1000.0e+12];   % Young's moduli
    case 12
        poissons = [0.0 0.0 0.0];           % poissons ratios
        youngs = [1000.0e+5 1000.0 1000.0];   % Young's moduli
    case 13
%       poissons = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ...
%             0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ...
%             0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ...
%             0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ...
%             0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ...
%             0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0];
%       poissons = [0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 ...
%             0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 ...
%             0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 ...
%             0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 ...
%             0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 ...
%             0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3];           % poisson's ratios
      poissons = [0.2962 0.3085 0.2986 0.2937 0.3081 0.3096 0.2988 ...
        0.2922 0.2952 0.2982 0.3019 0.2952 0.3021 0.3042 0.2944 0.2923 ...
        0.2959 0.2964 0.2985 0.3002 0.2917 0.2952 0.3060 0.2906 0.3086 ...
        0.3046 0.2998 0.3016 0.2947 0.2992 0.3093 0.3009 0.3004 0.2946 ...
        0.2998 0.3025 0.3036 0.2979 0.2973 0.3098 0.2908 0.3077 0.3083 ...
        0.3059 0.2920 0.2952 0.2967 0.3036 0.2927 0.3044];%0.29 + 0.02 * rand(1,50);
%       youngs = [210000 210000 210000 210000 210000 210000 210000 210000 ...
%         210000 210000 210000 210000 210000 210000 210000 210000 210000 ...
%         210000 210000 210000 210000 210000 210000 210000 ...
%         210000 210000 210000 210000 210000 210000 210000 210000 ...
%         210000 210000 210000 210000 210000 210000 210000 210000 210000 ...
%         210000 210000 210000 210000 210000 210000 210000 ...
%         210000 210000 210000 210000 210000 210000 210000 210000 ...
%         210000 210000 210000 210000 210000 210000 210000 210000 210000 ...
%         210000 210000 210000 210000 210000 210000 210000];
      youngs = 1.0e+5 * [2.0214 2.1308 2.0988 2.1558 2.1430 2.1807 ...
        2.1782 2.0668 2.1397 2.0396 2.0061 2.1488 2.1000 2.0960 2.1809 ...
        2.1220 2.1235 2.1719 2.1611 2.1153 2.0366 2.0480 2.1773 2.0057 ...
        2.0980 2.0336 2.1957 2.1425 2.1001 2.0942 2.0119 2.1364 2.0085 ...
        2.0143 2.1043 2.0193 2.1636 2.1635 2.1445 2.0300 2.1319 2.1037 ...
        2.1946 2.1298 2.1601 2.0908 2.0865 2.1651 2.0167 2.0266];%;200000 + 20000 * rand(1,50);
    otherwise
        error('MATLAB:preprocess:MaterialProperties',...
            'Unvalid Material ID "IFMatSet". Check ID or add additional case to switch-case-structure in "MaterialProperties.m".');
end;

MatProp.poissons = poissons;
MatProp.youngs = youngs;

clear poissons youngs;

end