function [stresse] = post_process(node,x,y,e,disp)
%2D quad element stress computation routine

global GRAININFO_ARR SUBELEMENT_GRAIN_MAP
xe = [];
ye = [];
dispj = [];
stressj = 0;

% grain = SUBELEMENT_GRAIN_MAP(e);

grain = 1;

% plane stress D matrix
%if grain == -1
%    D = [0 0 0;0 0 0;0 0 0];
%else
    young = GRAININFO_ARR(grain).youngs;
    pr = GRAININFO_ARR(grain).poisson;
    fac = young/(1 - (pr)^2);
    D = fac*[1.0, pr, 0;
            pr, 1.0, 0.0;
            0, 0, (1.-pr)/2 ];
%end

      
% get coordinates of element nodes 
for j=1:3
    je = node(j,e); xe(j) = x(je); ye(j) = y(je);
end
     
% element displacement vector
for j=1:3
    m1 = node(j,e)*2 - 1;
    m2 = node(j,e)*2;
    dispj(j*2-1) = disp(m1);
    dispj(j*2) = disp(m2);
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
        
% Area of the element
% Area = det([[1 1 1]' xe' ye'])/2;

%stress at centroid
stressj = stressj + (D*BJ*dispj')';

stresse(1,2:3) = xy;
stresse(1,4:6) = stressj;
stresse(1,1) = e;

