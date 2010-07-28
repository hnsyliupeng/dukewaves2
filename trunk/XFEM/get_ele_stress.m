% get_ele_stress.m
%
% CALL: get_ele_stress(elenodes,xcoords,ycoords, ...
%                        dis,disp_vec,BJ,cut,maxngrains)
%
% a routine for computing strains 'straine' and stresses
% 'stresse' elementwise, that computes the these quantities with respect
% to the subelements, too.
%
% Input paramters:
%   elenodes        global node IDs of current element
%   xcoords         x-coordinates of element's nodes
%   ycoords         y-coordinates of element's nodes
%   dis             global displacement vector (re-assembled)
%   disp_vec        global solution vector for element's DOFs
%   BJ              B-matrix for current element
%   cut             entry in 'cutlist' of current element
%   maxngrains      number of grains in discretization (Minimum: 3)
%
% Returned variables:
%   straine         strain in element 'e'
%   stresse         stress in element 'e'
%
% Structure of 'stresse':
%   Dimension: 1x6xmaxngrains
%   Index 1:    'stresse' for one element
%   Index 2:    column 1    global element ID
%               column 2    x-coordinate of element centroid
%               column 3    y-coordinate of element centroid
%               column 4    xx-stress at element centroid
%               column 5    yy-stress at element centroid
%               column 6    xy-stress at element centroid
%   Index 3:    ID of grain, to which these values belong to (maxngrains =
%               maximum number of grains)
%
% 'stresse' will be assembled to global stress matrix 'stress' in
% 'main_xfem_m'. There, it will have the same structure, but the first
% index will go up to the number of elements 'numele'.
%
% 'straine' is structured and used in a similar way.
%

% Author: Matthias Mayr (07/2010)

function [straine,stresse] = get_ele_stress(elenodes,xcoords,ycoords, ...
  dis,disp_vec,BJ,cut,maxngrains,e,id_dof)

% get access to global variables
global GRAININFO_ARR SUBELEMENT_GRAIN_MAP
global SUBELEM_INFO

% initialize some variables
dispj = [0 0 0 0 0 0];
strainj = [0 0 0];
stressj = [0 0 0];
basestrainj = [0 0 0];
enrstrainj = [0 0 0];
straine = zeros(1,6,maxngrains+1);
stresse = zeros(1,6,maxngrains);

% compute coordinates of centroid
centr = ones(1,3)/3;
xcen = centr*xcoords';
ycen = centr*ycoords';
xy = [xcen,ycen];

if cut == 0      % element is not cut by an interface
    grain = SUBELEMENT_GRAIN_MAP(e);
    young = GRAININFO_ARR(grain).youngs;
    pr = GRAININFO_ARR(grain).poisson;
    fac = young/(1 - (pr)^2);
    D = fac*[1.0, pr, 0;
          pr, 1.0, 0.0;
          0, 0, (1.-pr)/2 ];

    % element displacement vector
    for j=1:3
        m1 = elenodes(j)*2 - 1;   % x-direction
        m2 = elenodes(j)*2;       % y-direction
        dispj(j*2-1) = dis(m1);   % x-direction
        dispj(j*2) = dis(m2);     % y-direction
    end
    
    % strain and stress at centroid
    strainj = strainj + (BJ*dispj')';
    stressj = stressj + (D*BJ*dispj')';
    straine(1,2:3,1) = xy;
    straine(1,4:6,1) = strainj;
    straine(1,1,1) = e;
    stresse(1,2:3,grain) = xy;
    stresse(1,4:6,grain) = stressj;
    stresse(1,1,grain) = e;
else                    % element is cut by an interface
    % Calculate gradient of base degrees of freedom
    % Base element displacement vector
    for j=1:3
        m = elenodes(j);
        dispj(j*2-1) = disp_vec(m,1);
        dispj(j*2) = disp_vec(m,2);
    end

    %stress at centroid
    basestrainj = basestrainj + (BJ*dispj')';
    straine(1,2:3,1) = xy;
    straine(1,4:6,1) = basestrainj;
    straine(1,1,1) = e;
    dispj = [0 0 0 0 0 0];
    % Enriched element displacement vector
    for g = 1:maxngrains
        for j = 1:3
            m = elenodes(j);
            if find(id_dof(j) == g)
                find(id_dof(j) == g);
                dispj(j*2-1:j*2) = disp_vec(m,find(id_dof(j) == g));
            end
        end

        %stress at centroid
        enrstrainj = enrstrainj + (BJ*dispj')';
        straine(1,2:3,g+1) = xy;
        straine(1,4:6,g+1) = enrstrainj;
        straine(1,1,g+1) = e;
        dispj = [0 0 0 0 0 0];
        enrstrainj = 0;
    end
    D = zeros(3,3,maxngrains);

    % get all subelements of element 'e'
    subeleids = SUBELEM_INFO(1,e).kids;
    % get all grains, that overlap with element 'e'
    grains_of_ele = SUBELEMENT_GRAIN_MAP(subeleids);
    % loop over all grains
    for j = 1:maxngrains
        if any(grains_of_ele == j)  % only, if element 'e' overlaps with grain 'j'
            % get material parameters / build constitutive tensor
            young(j) = GRAININFO_ARR(j).youngs;
            pr(j) = GRAININFO_ARR(j).poisson;
            fac(j) = young(j)/((1) - (pr(j))^2);
            D(1:3,1:3,j) = fac(j)*[1.0, pr(j), 0;
                  pr(j), 1.0, 0.0;
                  0, 0, (1.-pr(j))/2 ];
            % compute stresses
            stresse(1,2:3,j) = xy;
            stresse(1,4:6,j) = D(:,:,j)*straine(1,4:6,1)' ...
                             + D(:,:,j)*straine(1,4:6,j+1)';
            stresse(1,1,j) = e;

        else    % element 'e' doesn't overlap with grain 'j'. So there
                % can't be stresses in element 'e', associated with that
                % grain 'j'. Just set them to zero.
            stresse(1,1:6,j) = [0 0 0 0 0 0];
        end;
    end
end

