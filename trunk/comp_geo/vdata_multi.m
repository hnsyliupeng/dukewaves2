% vdata_multi.m
% 
% This file contains datasets for the matrix 'p', which is needed for the
% voroni-tesselation. The variable "p" is 2 x (number of grains).  p(i,j) 
% gives the ith component of the center of grain j. The choice of the
% dataset, that will be used, is done by setting the parameter 'IFdatasetp'
% in the inpute file routine.
%
% To fill in an additional dataset in this file, you have to enlarge the
% switch-case-structure.

switch IFdatasetp
    case 1
        p = [0 0; 9 0; 4.5 9];      % 2 grains
    case 2
        p = [0 0; 13 1; 7.4 15];    % 2 grains, unsymmetric  
    case 3
        p = [0 -0.2; 3 0.2; 1.5 5]; % Debug problem
    case 4
        p = [0 0; 10 0; 5 9];       % 2 grains for beambending example
    case 5
        %Multi 6
        p = [ 0.0  2.1;
              0.1  0.2;
              0.0 -2.0;
              2.0  1.1;
              2.0 -1.0;
              3.9  2.0;
              4.0 -0.1;
              4.0 -2.3;
              6.1  1.0;
              6.0 -1.1;
              8.0  2.1;
              8.0  0.2;
              7.9 -2.0;
              9.9  1.1;
             10.2 -1.1;
             12.2  2.0;
             12.0 -0.1;
             11.9 -2.1;
             14.1  1.0;
             13.8 -0.9;
             16.1  2.2;
             16.0  0.2;
             15.9 -2.0];
    case 6
        %Multi1
        p =     [0    2.0000;
                 0    0.2000;
                 0   -2.0000;
            2.0000    1.1000;
            2.0000   -1.0000;
            4.0000    2.0000;
            4.0000   -0.1000;
            4.0000   -2.0000;
            6.0000    1.0500;
            6.0000   -1.1000;
            8.2000    2.0000;
            8.0000    0.1000;
            8.0000   -2.0000;
           10.0000    1.0000;
           10.0000   -1.1000;
           12.0000    2.0000;
           12.0000   -0.1000;
           12.0000   -2.0000;
           14.0000    1.0000;
           14.0000   -1.1000;
           16.0000    2.0000;
           16.0000    0.2000;
           16.0000   -2.0000];
    case 7
        % Multi5
        p =     [0    2.1000;
            0.1000    0.2000;
                 0   -2.0000;
            2.0000    1.1000;
            2.0000   -1.0000;
            3.9000    2.0000;
            4.0000   -0.1000;
            4.0000   -2.3000;
            6.1000    1.0000;
            6.0000   -1.1000;
            8.0000    2.1000;
            8.0000    0.2000;
            7.9000   -2.0000;
           10.1000    1.1000;
           10.0000   -1.2000;
           12.2000    2.0000;
           12.0000   -0.1000;
           11.9000   -2.1000;
           14.1000    1.0000;
           13.8000   -0.9000;
           16.1000    2.2000;
           16.0000    0.2000;
           16.0000   -2.0000];
    case 8
        %Multi4
        p =     [0    2.0000;
            0.1000    0.2000;
                 0   -2.0000;
            2.0000    1.1000;
            2.0000   -1.0000;
            3.9000    2.0000;
            4.0000   -0.1000;
            4.0000   -2.3000;
            6.1000    1.0000;
            6.0000   -1.1000;
            8.0000    2.1000;
            8.0000    0.2000;
            7.9000   -2.0000;
           10.1000    1.0000;
           10.0000   -1.1000;
           12.2000    2.0000;
           12.0000   -0.1000;
           12.0000   -2.1000;
           14.1000    1.0000;
           14.0000   -0.9000;
           16.1000    2.2000;
           16.0000    0.2000;
           16.0000   -2.0000];
    case 9
        % Multi3
        p =     [0    2.0000;
            0.1000    0.2000;
                 0   -2.0000;
            2.0000    1.1000;
            2.0000   -1.0000;
            3.9000    2.0000;
            4.0000   -0.1000;
            4.0000   -2.3000;
            6.1000    1.0000;
            6.0000   -1.1000;
            8.0000    2.1000;
            8.0000    0.2000;
            7.9000   -2.0000;
           10.1000    1.0000;
           10.0000   -1.1000;
           12.0000    2.0000;
           12.2000   -0.1000;
           12.0000   -2.1000;
           14.1000    1.0000;
           14.0000   -1.1000;
           16.1000    2.2000;
           16.0000    0.2000;
           16.1000   -2.0000];
    case 10
        % Multi2
        p =     [0    2.0000;
                 0    0.2000;
                 0   -2.0000;
            2.4000    1.1000;
            2.0000   -1.1000;
            3.9000    2.0000;
            4.0000   -0.1000;
            4.0000   -2.3000;
            6.1000    1.0000;
            5.9000   -1.1000;
            8.0000    2.1000;
            8.0000    0.2000;
            7.7000   -2.1000;
           10.0000    1.0000;
           10.0000   -1.1000;
           12.0000    2.0000;
           12.2000   -0.1000;
           12.0000   -2.0000;
           14.0000    1.0000;
           14.0000   -0.9000;
           16.0000    2.0000;
           16.0000    0.0000;
           16.0000   -2.0000];
    case 11
        p = [0 -2.6; 10 -2.6; 5 6.4];       % 3 grains
    case 12
        p = [-4 -3.1; 3.5 -3.1; 6 5.9];       % 3 grains
    case 13
        p = [9.1 0; 19.1 0; 14.1 9];       % 2 grains (vertical interface)
    case 14
%         p = [-1.0 0.0; 20.0 3.0; -1.0 6.0];    % 3 grains for
                                                % Hertzian_Contact_1
        p = [-1.0 0.0000001; 20.0 3.0000001; -1.0 6.0000001];
    case 15
        p = [-3 -2.6; 7 -2.6; 2 6.4];       % 3 grains
    case 16
        p = [10.5 1.5;12.5 -0.5;12.5 3.5];       % 3 grains
    case 17
        p = [-1.0 0.05; 20.0 3.05; -1.0 6.05];    % 3 grains for
                                                % Hertzian_Contact_2
    case 18
        p = [12.0 -2.0; 41.0 0.0; 12 2.0];    % 3 grains for
                                                % frictionless_sliding_anal
                                                % yt_1
    case 19
        p = [100 0; 110 0; 105 9]; % interfaces far away from domain --> only one grain
    case 20
        p = [12.0 -2.0; 41.0 0.0; 12 2.00001];    % 3 grains for
                                                % frictionless_sliding_anal
                                                % yt_1
    case 21
        p = [-1.0 -2.0; 100.0 0.0; -1 2.00001];    % 3 grains for
                                                % frictionless_sliding_anal
                                                % yt_3
    case 22
        p = [-1.0 -2.0; 100.0 0.0; -1 2.0];    % 3 grains for
                                                % frictionless_sliding_anal
                                                % yt_4
    case 23
        p = [-1.0 0.0; 6.0 0.5; -1 1.0];    % 3 grains for square
    case 24
        p = [2.0 -1.0; 3.0 -1; 4 5.0];    % 3 grains for beam_TM2
    case 25
        p = [0 0; 16 0; 8 9];       % 2 grains for frictionless_sliding_analyt_6
    case 26
        p = [0 3; 3*sin(pi/4) 3*cos(pi/4); 10 20];       % 2 grains for frictionless_sliding_analyt_7
    case 27
        p = [-23 0; 18 15; 18 -15];       % 2 grains for frictionless_sliding_analyt_8
    case 28
        p = [6 -2; 8 3; 10 -2];   
    case 29
        p = [0.2 -0.55; -0.8 1.05; 0.2 2.65];   % 3 grains for frictionless_sliding_analyt_2
                                                % unstructured mesh
    case 30
        p = [0.2 -1.55; -0.8 0.05; 0.2 1.65];   % 3 grains for frictionless_sliding_analyt_2
                                                % structured mesh
    otherwise
        error('MATLAB:comp_geo:vdata_multi',...
            'Unvalid ID for "IFdatasetp" !!!');
end;