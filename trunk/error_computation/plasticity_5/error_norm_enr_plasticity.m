% Calucates the L2 norm of the error in a problem with tri elements.  
% Jessica Sanders

disp('Compute error of displacement field for plasticity ...');

% use the following code lines to create a mat-file, which contains all
% data that are needed to interpolate the reference solution
%{
node_ref = node;
x_ref = x_orig;
y_ref = y_orig;
dis_ref = dis;
NODEINFO_ARR_ref = NODEINFO_ARR;
save C:\a_Daten\TUM\DA_Ausland\Matlab\numerical_results\plasticity_5_20100708\reference_solution\refsolu.mat ...
  node_ref x_ref y_ref dis_ref NODEINFO_ARR_ref;
save C:\a_Daten\TUM\DA_Ausland\Matlab\SVN_Code\error_computation\plasticity_5\refsolu.mat ...
  node_ref x_ref y_ref dis_ref NODEINFO_ARR_ref;
%}

% load reference solution
load refsolu.mat;

%Initialize error
total_errorx = 0;
total_errory = 0;

total_apprx = 0;
total_appry = 0;

error_in_elements = zeros(2,numele);  % stores error for each element to 
                                      % plot the error distribution

% Define 12-point Gauss quadrature
    
gr =    [0.873821971016996 0.063089014491502 0.063089014491502];
gs =    [0.063089014491502 0.063089014491502 0.873821971016996];
wg =    [0.050844906370207 0.050844906370207 0.050844906370207];
    
gr = [gr 0.501426509658179 0.249286745170910 0.249286745170910];
gs = [gs 0.249286745170910 0.249286745170910 0.501426509658179];
wg = [wg 0.116786275726379 0.116786275726379 0.116786275726379];
    
gr = [gr 0.636502499121399 0.310352451033785 0.053145049844816 0.636502499121399 0.310352451033785 0.053145049844816];
gs = [gs 0.310352451033785 0.053145049844816 0.636502499121399 0.053145049844816 0.636502499121399 0.310352451033785];
wg = [wg 0.082851075618374 0.082851075618374 0.082851075618374 0.082851075618374 0.082851075618374 0.082851075618374];


for j = 1:numele
    
    if (cutlist(j) == 0) % If element is uncut

        % initialize element error
        element_errorx = 0;
        element_errory = 0;
    
        approx_solnx = 0;
        approx_solny = 0;

        % get nodes:
        nodes = node(:,j);

        % get coordinates:
        xe = x(nodes);
        ye = y(nodes);
        
        % Area
        Area = det([[1 1 1]' xe' ye'])/2;
        
        % Shape function derivatives
        Ndr(1) = 1;
        Nds(1) = 0;
        Ndr(2) = 0;
        Nds(2) = 1;
        Ndr(3) = -1;
        Nds(3) = -1;
    
        % compute derivatives of x and y wrt psi and eta
        xdr = Ndr*xe'; ydr = Ndr*ye'; xds = Nds*xe';  yds = Nds*ye';
        jcob = xdr*yds - xds*ydr;

        % get nodal displacements
        xd = dis(2*nodes-1);
        yd = dis(2*nodes);

        % loop over Gauss points
        for i = 1:12
    
            r = gr(i);
            s = gs(i);
    
            % Shape functions
            N1 = r;
            N2 = s;
            N3 = 1-r-s; 

            % The displacements as a function of r and s
            dx_feta = N1*xd(1) + N2*xd(2) + N3*xd(3);
            dy_feta = N1*yd(1) + N2*yd(2) + N3*yd(3);
            
            % subtract the prescibed displacement
            dx_feta = dx_feta - 0.005;
    
            % x and y as a function of eta (for the analytical solution)
            x_feta = N1*xe(1) + N2*xe(2) + N3*xe(3);%- 0.005;
            y_feta = N1*ye(1) + N2*ye(2) + N3*ye(3) - 0.5;
            
            if y_feta > - 2   % check, in which grain the point resides
              % point (x_feta|y_feta) belongs to elastic block
              
              % Get analytical solution (calculate it here and not in a
              % function to accelerate the computation)

              % find node, which is the neares to current gauss point
              % Difference-vecotr
              diff_x = x_ref - x_feta;
              diff_y = y_ref - y_feta;

              % norm of difference
              norm_vec = sqrt(diff_x.^2 + diff_y.^2);

              % find node 'nodeID' with minimal distance
              [val nodeID] = min(norm_vec);

              % get elements connected to node 'nodeID'  
              elevec_ref = NODEINFO_ARR_ref(nodeID).elements;

              % initialize 'eleID'
              eleID_ref = 0;

              % Find element, which contains the current gauss point
              % loop over all elements, adjacent to 'nodeID'            
              for aa = elevec_ref
                % get nodes of current elements
                elenodes_ref = node_ref(:,aa);

                % coordinates of the element's nodes
                v = [x_ref(elenodes_ref); y_ref(elenodes_ref)];

                % point, which has to be checked
                p = [x_feta y_feta];

                % check, if current element 'aa' containts the gauss point
                % 'p'
                if polygon_contains_point_2d (3, v, p );
                  % 'aa' contains gauss point 'p'
                  eleID_ref = aa;
                end;
              end; 

              % get node IDs of element, that contains 'p'
              elenodes_ref = node_ref(:,eleID_ref);

              % get nodal coordinates
              xcoords_ref = x_ref(elenodes_ref);
              ycoords_ref = y_ref(elenodes_ref);

              % Interpolate reference solution at gauss point 'p' in element
              % 'eleID'
              
              % get area of reference element
              Area_ref = det([[1 1 1]' xcoords_ref' ycoords_ref']) / 2;
              
              % evaluate shape functions at gauss point 'p'
%               N_ref = zeros(2,6);
              N_ref2 = zeros(1,3);
              for b = 1:3     % Evaluate shape functions
                % load node coordinates
                xes = xcoords_ref;
                yes = ycoords_ref;

                % Get coordinates of area opposite node of concern
                xes(b) = p(1); 
                yes(b) = p(2);
                Larea_ref = det([[1 1 1]' xes' yes']) / 2;

%                 % Evaluate shape function for node 'b' (only in base DOFs)
%                 N_ref(1,2*b-1)  = N_ref(1,2*b-1)  + Larea_ref / Area_ref;
%                 N_ref(2,2*b)    = N_ref(2,2*b)    + Larea_ref / Area_ref;
                
                N_ref2(b) = N_ref2(b) + Larea_ref / Area_ref;
              end;

%               % get nodal displacements of 'eleID'
%               nodaldis_ref = dis_ref([2*elenodes_ref(1) - 1:2*elenodes_ref(1) ...
%                                      (2*elenodes_ref(2) - 1):2*elenodes_ref(2) ...
%                                       2*elenodes_ref(3) - 1:2*elenodes_ref(3)]);

              nodaldis_ref_x = dis_ref(2*elenodes_ref - 1);
              nodaldis_ref_y = dis_ref(2*elenodes_ref);
                                    
              % compute interpolation of reference displacement at gauss point
%               anal = N_ref * nodaldis_ref;
              anal(1) = N_ref2 * nodaldis_ref_x;
              anal(2) = N_ref2 * nodaldis_ref_y;

              % Difference between analytical and numerical solutions - error!
              ex = abs(dx_feta) - abs(anal(1));
              ey = abs(dy_feta) - abs(anal(2));

              % square the error
              ex2 = ex^2;
              ey2 = ey^2;

              % Assemble the squared error over the element
              element_errorx = element_errorx + ex2*wg(i)*jcob;
              element_errory = element_errory + ey2*wg(i)*jcob;

              % Assemble the approximate soln norm over the element
  %             approx_solnx = approx_solnx + (dx_feta)^2*wg(i)*jcob;
  %             approx_solny = approx_solny + (dy_feta)^2*wg(i)*jcob;
              approx_solnx = approx_solnx + (anal(1))^2*wg(i)*jcob;
              approx_solny = approx_solny + (anal(2))^2*wg(i)*jcob;
            else
              % point (x_feta|y_feta) belongs to rigid block
              anal(1) = 0;  % x-dispylacement of rigid block is zero
              anal(2) = 0;  % x-displacement of rigid block is zero
            end;
        end
    
        total_errorx = total_errorx + 0.5*element_errorx;
        total_errory = total_errory + 0.5*element_errory;
    
        total_apprx = total_apprx + 0.5*approx_solnx;
        total_appry = total_appry + 0.5*approx_solny;   
        
    elseif (cutlist(j) ~= 0) % If element is cut
        
        
        for m = 1:SUBELEM_INFO(j).no_kids % loop through subelements
        
            subele = SUBELEM_INFO(j).kids(m);
        
%             if (SUBELEMENT_GRAIN_MAP(subele) == 2) % If grain is not
%             enriched
            if (SUBELEMENT_GRAIN_MAP(subele) == 3) % If grain is not enriched
            
                % initialize element error
                element_errorx = 0;
                element_errory = 0;
    
                approx_solnx = 0;
                approx_solny = 0;

                % get nodes:
                node_sub = CONN(:,subele);

                % get coordinates:
                xe = X(node_sub);
                ye = Y(node_sub);
                
                % Shape function derivatives
                Ndr(1) = 1;
                Nds(1) = 0;
                Ndr(2) = 0;
                Nds(2) = 1;
                Ndr(3) = -1;
                Nds(3) = -1;
    
                % compute derivatives of x and y wrt psi and eta
                xdr = Ndr*xe'; ydr = Ndr*ye'; xds = Nds*xe';  yds = Nds*ye';
                jcob = xdr*yds - xds*ydr;
    
                Area = det([[1 1 1]' xe' ye'])/2;
                
                % get nodes of parent element
                nodes = node(:,j);

                % get nodal displacements at parent nodes
                xd = dis(2*nodes-1);
                yd = dis(2*nodes);    

                % loop over Gauss points
                for i = 1:12
    
                    r = gr(i);
                    s = gs(i);
    
                    % Shape functions
                    N1 = r;
                    N2 = s;
                    N3 = 1-r-s; 
        
%                     % Shape function derivatives
%                     Ndr(1) = 1;
%                     Nds(1) = 0;
%                     Ndr(2) = 0;
%                     Nds(2) = 1;
%                     Ndr(3) = -1;
%                     Nds(3) = -1;
%     
%                     % compute derivatives of x and y wrt psi and eta
%                     xdr = Ndr*xe'; ydr = Ndr*ye'; xds = Nds*xe';  yds = Nds*ye';
%                     jcob = xdr*yds - xds*ydr;
    
%                     Area = det([[1 1 1]' xe' ye'])/2;

                    % The displacements as a function of r and s
                    dx_feta = N1*xd(1) + N2*xd(2) + N3*xd(3);
                    dy_feta = N1*yd(1) + N2*yd(2) + N3*yd(3);
                    
                    dx_feta = dx_feta - 0.005;
    
                    % x and y as a function of eta (for the analytical solution)
                    x_feta = N1*xe(1) + N2*xe(2) + N3*xe(3);
                    y_feta = N1*ye(1) + N2*ye(2) + N3*ye(3) - 0.5;
    
                    if y_feta > - 2   % check, in which grain the point resides
                      % point (x_feta|y_feta) belongs to elastic block

                      % Get analytical solution (calculate it here and not in a
                      % function to accelerate the computation)

                      % find node, which is the neares to current gauss point
                      % Difference-vecotr
                      diff_x = x_ref - x_feta;
                      diff_y = y_ref - y_feta;

                      % norm of difference
                      norm_vec = sqrt(diff_x.^2 + diff_y.^2);

                      % find node 'nodeID' with minimal distance
                      [val nodeID] = min(norm_vec);

                      % get elements connected to node 'nodeID'  
                      elevec_ref = NODEINFO_ARR_ref(nodeID).elements;

                      % initialize 'eleID'
                      eleID_ref = 0;

                      % Find element, which contains the current gauss point
                      % loop over all elements, adjacent to 'nodeID'            
                      for aa = elevec_ref
                        % get nodes of current elements
                        elenodes_ref = node_ref(:,aa);

                        % coordinates of the element's nodes
                        v = [x_ref(elenodes_ref); y_ref(elenodes_ref)];

                        % point, which has to be checked
                        p = [x_feta y_feta];

                        % check, if current element 'aa' containts the gauss point
                        % 'p'
                        if polygon_contains_point_2d (3, v, p );
                          % 'aa' contains gauss point 'p'
                          eleID_ref = aa;
                        end;
                      end; 

                      % get node IDs of element, that contains 'p'
                      elenodes_ref = node_ref(:,eleID_ref);

                      % get nodal coordinates
                      xcoords_ref = x_ref(elenodes_ref);
                      ycoords_ref = y_ref(elenodes_ref);

                      % Interpolate reference solution at gauss point 'p' in element
                      % 'eleID'

                      % get area of reference element
                      Area_ref = det([[1 1 1]' xcoords_ref' ycoords_ref']) / 2;

                      % evaluate shape functions at gauss point 'p'
        %               N_ref = zeros(2,6);
                      N_ref2 = zeros(1,3);
                      for b = 1:3     % Evaluate shape functions
                        % load node coordinates
                        xes = xcoords_ref;
                        yes = ycoords_ref;

                        % Get coordinates of area opposite node of concern
                        xes(b) = p(1); 
                        yes(b) = p(2);
                        Larea_ref = det([[1 1 1]' xes' yes']) / 2;

        %                 % Evaluate shape function for node 'b' (only in base DOFs)
        %                 N_ref(1,2*b-1)  = N_ref(1,2*b-1)  + Larea_ref / Area_ref;
        %                 N_ref(2,2*b)    = N_ref(2,2*b)    + Larea_ref / Area_ref;

                        N_ref2(b) = N_ref2(b) + Larea_ref / Area_ref;
                      end;

        %               % get nodal displacements of 'eleID'
        %               nodaldis_ref = dis_ref([2*elenodes_ref(1) - 1:2*elenodes_ref(1) ...
        %                                      (2*elenodes_ref(2) - 1):2*elenodes_ref(2) ...
        %                                       2*elenodes_ref(3) - 1:2*elenodes_ref(3)]);

                      nodaldis_ref_x = dis_ref(2*elenodes_ref - 1);
                      nodaldis_ref_y = dis_ref(2*elenodes_ref);

                      % compute interpolation of reference displacement at gauss point
        %               anal = N_ref * nodaldis_ref;
                      anal(1) = N_ref2 * nodaldis_ref_x;
                      anal(2) = N_ref2 * nodaldis_ref_y;

                      % Difference between analytical and numerical solutions - error!
                      ex = abs(dx_feta) - abs(anal(1));
                      ey = abs(dy_feta) - abs(anal(2));

                      % square the error
                      ex2 = ex^2;
                      ey2 = ey^2;

                      % Assemble the squared error over the element
                      element_errorx = element_errorx + ex2*wg(i)*jcob;
                      element_errory = element_errory + ey2*wg(i)*jcob;

                      % Assemble the approximate soln norm over the element
          %             approx_solnx = approx_solnx + (dx_feta)^2*wg(i)*jcob;
          %             approx_solny = approx_solny + (dy_feta)^2*wg(i)*jcob;
                      approx_solnx = approx_solnx + (anal(1))^2*wg(i)*jcob;
                      approx_solny = approx_solny + (anal(2))^2*wg(i)*jcob;
                    else
                      % point (x_feta|y_feta) belongs to rigid block
                      anal(1) = 0;  % x-dispylacement of rigid block is zero
                      anal(2) = 0;  % x-displacement of rigid block is zero
                    end;
    
                  % Assemble the approximate soln norm over the element
                  approx_solnx = approx_solnx + (dx_feta)^2*wg(i)*jcob;
                  approx_solny = approx_solny + (dy_feta)^2*wg(i)*jcob;
        
                end

                total_errorx = total_errorx + 0.5*element_errorx;
                total_errory = total_errory + 0.5*element_errory;
    
                total_apprx = total_apprx + 0.5*approx_solnx;
                total_appry = total_appry + 0.5*approx_solny;

        
%             elseif (SUBELEMENT_GRAIN_MAP(subele) == 1)  % If grain is enriched
            elseif (SUBELEMENT_GRAIN_MAP(subele) == 1)  % If grain is enriched
            
                % initialize element error
                element_errorx = 0;
                element_errory = 0;
    
                approx_solnx = 0;
                approx_solny = 0;

                % get nodes of subelement
                node_sub = CONN(:,subele);

                % get coordinates of subelement nodes
                xe = X(node_sub);
                ye = Y(node_sub);
                
                % Shape function derivatives
                Ndr(1) = 1;
                Nds(1) = 0;
                Ndr(2) = 0;
                Nds(2) = 1;
                Ndr(3) = -1;
                Nds(3) = -1;
    
                % compute derivatives of x and y wrt psi and eta
                xdr = Ndr*xe'; ydr = Ndr*ye'; xds = Nds*xe';  yds = Nds*ye';
                jcob = xdr*yds - xds*ydr;
        
                Area = det([[1 1 1]' xe' ye'])/2;
                
                % get nodes of parent element
                nodes = node(:,j);
        
                % element displacement vector
                for b=1:3
                    % base degrees of freedom
                    b1 = id_eqns(node(b,j),1);
                    b2 = id_eqns(node(b,j),2); 
                    dispx(b) = fdisp(b1);
                    dispy(b) = fdisp(b2);
                end
        
                % element displacement vector
                for b=1:3
                    % extra degrees of freedom
                    b1 = id_eqns(node(b,j),3);
                    b2 = id_eqns(node(b,j),4); 
                    dispx(b) = dispx(b) + fdisp(b1);
                    dispy(b) = dispy(b) + fdisp(b2);
                end

                % get nodal displacements at parent nodes
                xd = dispx;%dis(2*nodes-1);
                yd = dispy;%dis(2*nodes);    
                
                % loop over Gauss points
                for i = 1:12
    
                    r = gr(i);
                    s = gs(i);
    
                    % Shape functions
                    N1 = r;
                    N2 = s;
                    N3 = 1-r-s; 

                    % The displacements as a function of r and s
                    dx_feta = N1*xd(1) + N2*xd(2) + N3*xd(3);
                    dy_feta = N1*yd(1) + N2*yd(2) + N3*yd(3);
                    
                    dx_feta = dx_feta - 0.005;%
    
                    % x and y as a function of eta (for the analytical solution)
                    x_feta = N1*xe(1) + N2*xe(2) + N3*xe(3);
                    y_feta = N1*ye(1) + N2*ye(2) + N3*ye(3) - 0.5;

                    if y_feta > - 2   % check, in which grain the point resides
                      % point (x_feta|y_feta) belongs to elastic block
                      
                      % Get analytical solution (calculate it here and not in a
                      % function to accelerate the computation)

                      % find node, which is the neares to current gauss point
                      % Difference-vecotr
                      diff_x = x_ref - x_feta;
                      diff_y = y_ref - y_feta;

                      % norm of difference
                      norm_vec = sqrt(diff_x.^2 + diff_y.^2);

                      % find node 'nodeID' with minimal distance
                      [val nodeID] = min(norm_vec);

                      % get elements connected to node 'nodeID'  
                      elevec_ref = NODEINFO_ARR_ref(nodeID).elements;

                      % initialize 'eleID'
                      eleID_ref = 0;

                      % Find element, which contains the current gauss point
                      % loop over all elements, adjacent to 'nodeID'            
                      for aa = elevec_ref
                        % get nodes of current elements
                        elenodes_ref = node_ref(:,aa);

                        % coordinates of the element's nodes
                        v = [x_ref(elenodes_ref); y_ref(elenodes_ref)];

                        % point, which has to be checked
                        p = [x_feta y_feta];

                        % check, if current element 'aa' containts the gauss point
                        % 'p'
                        if polygon_contains_point_2d (3, v, p );
                          % 'aa' contains gauss point 'p'
                          eleID_ref = aa;
                        end;
                      end; 

                      % get node IDs of element, that contains 'p'
                      elenodes_ref = node_ref(:,eleID_ref);

                      % get nodal coordinates
                      xcoords_ref = x_ref(elenodes_ref);
                      ycoords_ref = y_ref(elenodes_ref);

                      % Interpolate reference solution at gauss point 'p' in element
                      % 'eleID'

                      % get area of reference element
                      Area_ref = det([[1 1 1]' xcoords_ref' ycoords_ref']) / 2;

                      % evaluate shape functions at gauss point 'p'
        %               N_ref = zeros(2,6);
                      N_ref2 = zeros(1,3);
                      for b = 1:3     % Evaluate shape functions
                        % load node coordinates
                        xes = xcoords_ref;
                        yes = ycoords_ref;

                        % Get coordinates of area opposite node of concern
                        xes(b) = p(1); 
                        yes(b) = p(2);
                        Larea_ref = det([[1 1 1]' xes' yes']) / 2;

        %                 % Evaluate shape function for node 'b' (only in base DOFs)
        %                 N_ref(1,2*b-1)  = N_ref(1,2*b-1)  + Larea_ref / Area_ref;
        %                 N_ref(2,2*b)    = N_ref(2,2*b)    + Larea_ref / Area_ref;

                        N_ref2(b) = N_ref2(b) + Larea_ref / Area_ref;
                      end;

        %               % get nodal displacements of 'eleID'
        %               nodaldis_ref = dis_ref([2*elenodes_ref(1) - 1:2*elenodes_ref(1) ...
        %                                      (2*elenodes_ref(2) - 1):2*elenodes_ref(2) ...
        %                                       2*elenodes_ref(3) - 1:2*elenodes_ref(3)]);

                      nodaldis_ref_x = dis_ref(2*elenodes_ref - 1);
                      nodaldis_ref_y = dis_ref(2*elenodes_ref);

                      % compute interpolation of reference displacement at gauss point
        %               anal = N_ref * nodaldis_ref;
                      anal(1) = N_ref2 * nodaldis_ref_x;
                      anal(2) = N_ref2 * nodaldis_ref_y;

                      % Difference between analytical and numerical solutions - error!
                      ex = abs(dx_feta) - abs(anal(1));
                      ey = abs(dy_feta) - abs(anal(2));

                      % square the error
                      ex2 = ex^2;
                      ey2 = ey^2;

                      % Assemble the squared error over the element
                      element_errorx = element_errorx + ex2*wg(i)*jcob;
                      element_errory = element_errory + ey2*wg(i)*jcob;

                      % Assemble the approximate soln norm over the element
          %             approx_solnx = approx_solnx + (dx_feta)^2*wg(i)*jcob;
          %             approx_solny = approx_solny + (dy_feta)^2*wg(i)*jcob;
                      approx_solnx = approx_solnx + (anal(1))^2*wg(i)*jcob;
                      approx_solny = approx_solny + (anal(2))^2*wg(i)*jcob;
                    else
                      % point (x_feta|y_feta) belongs to rigid block
                      anal(1) = 0;  % x-dispylacement of rigid block is zero
                      anal(2) = 0;  % x-displacement of rigid block is zero
                    end;
        
                    % Assemble the approximate soln norm over the element
                    approx_solnx = approx_solnx + (dx_feta)^2*wg(i)*jcob;
                    approx_solny = approx_solny + (dy_feta)^2*wg(i)*jcob;
        
        
                end
    
                total_errorx = total_errorx + 0.5*element_errorx;
                total_errory = total_errory + 0.5*element_errory;
    
                total_apprx = total_apprx + 0.5*approx_solnx;
                total_appry = total_appry + 0.5*approx_solny;   
            end
        end
    end
    
    error_in_elements(1,j) = sqrt(element_errorx^2 + element_errory^2);
    error_in_elements(2,j) = sqrt(approx_solnx^2 + approx_solny^2);
end

% prepare a figure to plot the error distribution over the domain
figure(70)
hold on;
set(70,'Name','Error in displacement field (L2-Norm)');
minerror = min(error_in_elements(1,:));% ./ error_in_elements(2,:));
maxerror = max(error_in_elements(1,:));% ./ error_in_elements(2,:));
V = [minerror maxerror];
colormap(jet(100));
caxis(V);
colorbar;

for j=1:numele
  xcoords = x(node(:,j));
  ycoords = y(node(:,j));
  if error_in_elements(2,j) ~= 0
    ele_error = error_in_elements(1,j);% / error_in_elements(2,j);
  else
    ele_error = 0;
  end
  patch(xcoords,ycoords,ele_error,'EdgeColor','none');
end;

L2norm = sqrt((total_errorx + total_errory)/(total_apprx + total_appry));
disp(['L2-norm of displacement:   ' num2str(L2norm)]);
