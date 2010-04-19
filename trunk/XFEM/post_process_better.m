function [straine,stresse] = post_process_better(node,x,y,e,disp,disp_vec,id_dof,cutlist,maxngrains)
%2D quad element stress computation routine

global GRAININFO_ARR SUBELEMENT_GRAIN_MAP
global SUBELEM_INFO
xe = [];
ye = [];
dispj = [0 0 0 0 0 0];
strainj = [0 0 0];
stressj = [0 0 0];
basestrainj = [0 0 0];
enrstrainj = [0 0 0];
straine = zeros(1,6,maxngrains+1);
stresse = zeros(1,6,maxngrains);
% get coordinates of element nodes 
for j=1:3
    je = node(j,e); xe(j) = x(je); ye(j) = y(je);
end

% compute derivatives of shape functions in reference coordinates
NJr(1) = 1;
NJr(2) = 0;
NJr(3) = -1;
NJs(1) = 0;
NJs(2) = 1;
NJs(3) = -1;
% compute derivatives of x and y wrt psi and eta
xr = NJr*xe'; yr = NJr*ye'; xs = NJs*xe';  ys = NJs*ye';
Jinv = [ys, -yr; -xs, xr];
jcob = xr*ys - xs*yr;
% compute derivatives of shape functions in element coordinates
NJdrs = [NJr; NJs];
NJdxy = Jinv*NJdrs/jcob;
% assemble B matrix
BJ = zeros(3,6);
BJ(1,1:2:5) = NJdxy(1,1:3);  BJ(2,2:2:6) = NJdxy(2,1:3);
BJ(3,1:2:5) = NJdxy(2,1:3);  BJ(3,2:2:6) = NJdxy(1,1:3);

% compute coordinates of centroid
centr = ones(1,3)/3;
xcen = centr*xe';
ycen = centr*ye';
xy = [xcen,ycen];


if cutlist(e) == 0
    
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
        dispj(j*2-1) = disp(m1);
        dispj(j*2) = disp(m2);
    end    

    % Area of the element
    % Area = det([[1 1 1]' xe' ye'])/2;

    %stress at centroid
    strainj = strainj + (BJ*dispj')';
    stressj = stressj + (D*BJ*dispj')';

    straine(1,2:3,1) = xy;
    straine(1,4:6,1) = strainj;
    straine(1,1,1) = e;
    
    stresse(1,2:3,grain) = xy;
    stresse(1,4:6,grain) = stressj;
    stresse(1,1,grain) = e;
    
else
    
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
    
        dispj;
        
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
                m;
            dispj(j*2-1:j*2) = disp_vec(m,find(id_dof(m,:) == g));
            end
        end
            dispj;
            
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
    
    for j = 1:maxngrains

    young(j) = GRAININFO_ARR(j).youngs;
    pr(j) = GRAININFO_ARR(j).poisson;
    fac(j) = young(j)/((1) - (pr(j))^2);
    D(1:3,1:3,j) = fac(j)*[1.0, pr(j), 0;
          pr(j), 1.0, 0.0;
          0, 0, (1.-pr(j))/2 ];

    stresse(1,2:3,j) = xy;
    stresse(1,4:6,j) = D(:,:,j)*straine(1,4:6,1)' ...
                     + D(:,:,j)*straine(1,4:6,j+1)';
    stresse(1,1,j) = e;
    end
    
end




