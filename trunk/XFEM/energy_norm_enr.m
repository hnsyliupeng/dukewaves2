% Calucates the L2 norm of the error in a problem with tri elements.  
% Jessica Sanders

%Initialize error
total_error = 0;

total_appr = 0;

young = 1000;
pr = 0.3;
fac = young/(1 - (pr)^2);
D = fac*[1.0, pr, 0;
    pr, 1.0, 0.0;
    0, 0, (1.-pr)/2 ];


% compute derivatives of shape functions in reference coordinates
NJr(1) = 1;
NJr(2) = 0;
NJr(3) = -1;
NJs(1) = 0;
NJs(2) = 1;
NJs(3) = -1;

% Define 6-point Gauss quadrature
% gr =    [0.816847572980459 0.091576213509771 0.091576213509771];
% gs =    [0.091576213509771 0.091576213509771 0.816847572980459];
% wg =    [0.109951743655322 0.109951743655322 0.109951743655322];

% gr = [gr 0.108103018168070 0.445948490915965 0.445948490915965];
% gs = [gs 0.445948490915965 0.445948490915965 0.108103018168070];
% wg = [wg 0.223381589678011 0.223381589678011 0.223381589678011];

% Define 12-point Gauss quadrature

gr =    [0.873821971016996 0.063089014491502 0.063089014491502];
gs =    [0.063089014491502 0.063089014491502 0.873821971016996];
wg =    [0.050844906370207 0.050844906370207 0.050844906370207];

gr = [gr 0.501426509658179 0.249286745170910 0.249286745170910];
gs = [gs 0.249286745170910 0.249286745170910 0.501426509658179];
wg = [wg 0.116786275726379 0.116786275726379 0.116786275726379];

gr = [gr 0.636502499121399 0.310352451033785 0.053145049844816];
gs = [gs 0.310352451033785 0.053145049844816 0.636502499121399];
wg = [wg 0.082851075618374 0.082851075618374 0.082851075618374];

gr = [gr 0.636502499121399 0.310352451033785 0.053145049844816];
gs = [gs 0.053145049844816 0.636502499121399 0.310352451033785];
wg = [wg 0.082851075618374 0.082851075618374 0.082851075618374];

for j = 1:numele
    
    % initialize nodal displacements
    nodal_disp   = [];
    nodal_disp_e = [];
    
    % initialize element error
    element_error = 0;
    approx_soln = 0;
    
    % get nodes of parent element
    nodes = node(:,j);

    % get coordinates of parent element
    xe = x(nodes);
    ye = y(nodes);

    % get nodal displacements
    for i = 1:3
        nodal_disp(2*i-1) = dis(2*nodes(i) - 1);
        nodal_disp(2*i)   = dis(2*nodes(i));
    end

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
    % Area of the element
    Area = det([[1 1 1]' xe' ye'])/2;


    if (cutlist(j) == 0) % If element is uncut
        
        eps = BJ*nodal_disp';

        % loop over Gauss points
        for i = 1:12
    
            r = gr(i);
            s = gs(i);
    
            % Shape functions
            N1 = r;
            N2 = s;
            N3 = 1-r-s; 

            % x and y as a function of eta (for the analytical solution)
            x_feta = N1*xe(1) + N2*xe(2) + N3*xe(3);
            y_feta = N1*ye(1) + N2*ye(2) + N3*ye(3);
    
            % Get analytical solution
            anal = analytical_eng(x_feta, y_feta);
    
            % Difference between analytical and numerical solutions - error!
            e = eps - anal';
    
            % The quantity in the energy norm
            eng = e'*D*e;
        
            % The approx soln
            hsol = eps'*D*eps;
    
            % Assemble the squared error over the element
            element_error = element_error + eng*wg(i)*Area;
        
            % Assemble the approximate soln norm over the element
            approx_soln = approx_soln + hsol*wg(i)*Area;
                
        end
        
    elseif (cutlist(j) ~= 0) % If element is cut
        
        for m = 1:SUBELEM_INFO(j).no_kids % loop through subelements
        
            subele = SUBELEM_INFO(j).kids(m);
        
            if (SUBELEMENT_GRAIN_MAP(subele) == 2) % If grain is not enriched
                
                eps = BJ*nodal_disp';
                               
                % get nodes of subelement
                node_sub = CONN(:,subele);

                % get coordinates of 
                xes = X(node_sub);
                yes = X(node_sub);
                
                Area = det([[1 1 1]' xes' yes'])/2;
                
                % loop over Gauss points
                for i = 1:12

                    r = gr(i);
                    s = gs(i);

                    % Shape functions
                    N1 = r;
                    N2 = s;
                    N3 = 1-r-s;

                    % x and y as a function of eta (for the analytical solution)
                    x_feta = N1*xes(1) + N2*xes(2) + N3*xes(3);
                    y_feta = N1*yes(1) + N2*yes(2) + N3*yes(3);

                    % Get analytical solution
                    anal = analytical_eng(x_feta, y_feta);

                    % Difference between analytical and numerical solutions - error!
                    e = eps - anal';

                    % The quantity in the energy norm
                    eng = e'*D*e;

                    % The approx soln
                    hsol = eps'*D*eps;

                    % Assemble the squared error over the element
                    element_error = element_error + eng*wg(i)*Area;

                    % Assemble the approximate soln norm over the element
                    approx_soln = approx_soln + hsol*wg(i)*Area;

                end

            elseif (SUBELEMENT_GRAIN_MAP(subele) == 1)  % If grain is enriched 
                
                % get nodes of subelement
                node_sub = CONN(:,subele);

                % get coordinates of
                xes = X(node_sub);
                yes = X(node_sub);

                Area = det([[1 1 1]' xes' yes'])/2;
                     
                % element displacement vector
                for b=1:3
                    % base degrees of freedom
                    b1 = id_eqns(node(b,j),1);
                    b2 = id_eqns(node(b,j),2); 
                    nodal_disp(2*b-1) = fdisp(b1);
                    nodal_disp(2*b)   = fdisp(b2);
                end
        
                % element displacement vector
                for b=1:3
                    % extra degrees of freedom
                    b1 = id_eqns(node(b,j),3);
                    b2 = id_eqns(node(b,j),4); 
                    nodal_disp_e(2*b-1) = fdisp(b1);
                    nodal_disp_e(2*b)   = fdisp(b2);
                end
                
                eps = BJ*nodal_disp';
                
                eps_e = BJ*nodal_disp_e';
                
                eps = eps + eps_e;
                
                % loop over Gauss points
                for i = 1:12

                    r = gr(i);
                    s = gs(i);

                    % Shape functions
                    N1 = r;
                    N2 = s;
                    N3 = 1-r-s;

                    % x and y as a function of eta (for the analytical solution)
                    x_feta = N1*xes(1) + N2*xes(2) + N3*xes(3);
                    y_feta = N1*yes(1) + N2*yes(2) + N3*yes(3);

                    % Get analytical solution
                    anal = analytical_eng(x_feta, y_feta);

                    % Difference between analytical and numerical solutions - error!
                    e = eps - anal';

                    % The quantity in the energy norm
                    eng = e'*D*e;

                    % The approx soln
                    hsol = eps'*D*eps;

                    % Assemble the squared error over the element
                    element_error = element_error + eng*wg(i)*Area;

                    % Assemble the approximate soln norm over the element
                    approx_soln = approx_soln + hsol*wg(i)*Area;

                end
            end
        end
    end
    total_error = total_error + element_error;    
    total_appr = total_appr + approx_soln;    
end

energy_norm = sqrt(total_error)
