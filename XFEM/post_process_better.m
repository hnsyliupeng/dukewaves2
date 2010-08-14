% post_process_better.m
%
% a post processing routine for computing strains 'straine' and stresses
% 'stresse' elementwise, that computes the these quantities with respect
% to the subelements, too.
%
% Input paramters:
%   node            mapping between elements and their nodes
%   x               x-coordinates of all nodes
%   y               y-coordinates of all nodes
%   e               global element ID
%   dis             global displacement vector (re-assembled)
%   disp_vec
%   id_dof          stores, which DOFs are enriched once or twice
%   cutlist         lists, which DOFs are enriched
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
function [straine,stresse] = post_process_better(node,x,y,e,dis, ...
    disp_vec,id_dof,cutlist,maxngrains)
% get access to global variables
global GRAININFO_ARR SUBELEMENT_GRAIN_MAP
global SUBELEM_INFO
% initialize some variables
xe = zeros(1,3);
ye = zeros(1,3);
dispj = [0 0 0 0 0 0];
strainj = [0 0 0];
stressj = [0 0 0];
basestrainj = [0 0 0];
enrstrainj = [0 0 0];
straine = zeros(1,6,maxngrains+1);
stresse = zeros(1,6,maxngrains);
% get coordinates of nodes of element 'e'
for j=1:3
  je = node(j,e); 
  xe(j) = x(je); 
  ye(j) = y(je);
end
% compute derivatives of shape functions in reference coordinates
NJr(1) = 1;
NJr(2) = 0;
NJr(3) = -1;
NJs(1) = 0;
NJs(2) = 1;
NJs(3) = -1;
% compute derivatives of x and y wrt psi and eta
xr = NJr*xe'; 
yr = NJr*ye'; 
xs = NJs*xe';  
ys = NJs*ye';
Jinv = [ys, -yr; -xs, xr];
jcob = xr*ys - xs*yr;
% compute derivatives of shape functions in element coordinates
NJdrs = [NJr; NJs];
NJdxy = Jinv*NJdrs/jcob;
% assemble B matrix
BJ = zeros(3,6);
BJ(1,1:2:5) = NJdxy(1,1:3);  
BJ(2,2:2:6) = NJdxy(2,1:3);
BJ(3,1:2:5) = NJdxy(2,1:3);  
BJ(3,2:2:6) = NJdxy(1,1:3);
% compute coordinates of centroid
centr = ones(1,3)/3;
xcen = centr*xe';
ycen = centr*ye';
xy = [xcen,ycen];

if cutlist(e) == 0      % element is not cut by an interface
    grain = SUBELEMENT_GRAIN_MAP(e);
    young = GRAININFO_ARR(grain).youngs;
    pr = GRAININFO_ARR(grain).poisson;
    fac = young/(1 - (pr)^2);
    D = fac*[1.0, pr, 0;
          pr, 1.0, 0.0;
          0, 0, (1.-pr)/2 ];

    % element displacement vector
    for j=1:3
        m1 = node(j,e)*2 - 1;
        m2 = node(j,e)*2;
        dispj(j*2-1) = dis(m1);
        dispj(j*2) = dis(m2);
    end
    % Area of the element
    % Area = det([[1 1 1]' xe' ye'])/2;
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
%     for j=1:3
%         m = node(j,e);
%         dispj(j*2-1) = disp_vec(m,1);
%         dispj(j*2) = disp_vec(m,2);
%     end
    for j=1:3
        m = node(j,e);
        dispj(j*2-1) = disp_vec(m,1);
        dispj(j*2) = disp_vec(m,2);
    end
%     grain = SUBELEMENT_GRAIN_MAP(e);
%
%     young = GRAININFO_ARR(grain).youngs;
%     pr = GRAININFO_ARR(grain).poisson;
%     fac = young/(1 - (pr)^2);
%     D = fac*[1.0, pr, 0;
%           pr, 1.0, 0.0;
%           0, 0, (1.-pr)/2 ];

    %stress at centroid
    basestrainj = basestrainj + (BJ*dispj')';
    straine(1,2:3,1) = xy;
    straine(1,4:6,1) = basestrainj;
    straine(1,1,1) = e;
    dispj = [0 0 0 0 0 0];
    % Enriched element displacement vector
    for g = 1:maxngrains
        for j = 1:3
            m = node(j,e);
            if find(id_dof(m,:) == g)
                find(id_dof(m,:) == g);
                dispj(j*2-1:j*2) = disp_vec(m,find(id_dof(m,:) == g));
            end
        end

%         grain = g;
%
%         young = GRAININFO_ARR(grain).youngs;
%         pr = GRAININFO_ARR(grain).poisson;
%         fac = young/(1 - (pr)^2);
%         D = fac*[1.0, pr, 0;
%                 pr, 1.0, 0.0;
%                 0, 0, (1.-pr)/2 ];

        %stress at centroid
        enrstrainj = enrstrainj + (BJ*dispj')';
        straine(1,2:3,g+1) = xy;
        straine(1,4:6,g+1) = enrstrainj;
        straine(1,1,g+1) = e;
        dispj = [0 0 0 0 0 0];
        enrstrainj = 0;
    end
    D = zeros(3,3,maxngrains);
%     for j = 1:maxngrains      % from Jessica
%
%         young(j) = GRAININFO_ARR(j).youngs;
%         pr(j) = GRAININFO_ARR(j).poisson;
%         fac(j) = young(j)/((1) - (pr(j))^2);
%         D(1:3,1:3,j) = fac(j)*[1.0, pr(j), 0;
%               pr(j), 1.0, 0.0;
%               0, 0, (1.-pr(j))/2 ];
%
%         stresse(1,2:3,j) = xy;
%         stresse(1,4:6,j) = D(:,:,j)*straine(1,4:6,1)' ...
%                          + D(:,:,j)*straine(1,4:6,j+1)';
%         stresse(1,1,j) = e;
%
%         if stresse(1,5,3) < -7
%             save test.mat
%         end;
%     end
% here start Matthias' modifications
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
