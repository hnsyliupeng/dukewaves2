% displacement_error_reference.m
%
% Computes the L2-norm of the error in the displacement field. Since an
% analytical solution is not possible, a fine mesh solution is assumed to
% be the analytial solution. The convergence of this fine mesh solution has
% been shown at some characteristic nodes.

% Author: Matthias Mayr (07/2010)

disp('displacement_error_reference ...');

%% create dataset for reference solution
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
% ----------------------------------------------------------------------- %
%% Initialize
% load the reference solution
load refsolu.mat;

% set gauss points and weights for a 12-point-gauss-quadrature
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

% error
total_error = 0;

% exact
total_solu = 0;
% ----------------------------------------------------------------------- %
%% computation
% loop over all elements
for e = 1:numele
  ele_error = 0;
  ele_solu = 0;
  
  % get nodes of currente element 'e'
  elenodes  = node(:,e);
  
  % get x- and y-coordinates of nodes
  xcoords = x(elenodes);
  ycoords = y(elenodes);
  elecoords = [ xcoords(1) 
                ycoords(1)
                xcoords(2) 
                ycoords(2)
                xcoords(3) 
                ycoords(3)];
  
  % compute area of current element
  area = det([[1 1 1]' xcoords' ycoords']) / 2;
  
  % get element displacement vector
  eledis = [dis(2*elenodes(1) - 1)
            dis(2*elenodes(1))
            dis(2*elenodes(2) - 1)
            dis(2*elenodes(2))
            dis(2*elenodes(3) - 1)
            dis(2*elenodes(3))     ];
   
  % compute derivatives of shape functions in reference coordinates
  NJr(1) = 1;
  NJr(2) = 0;
  NJr(3) = -1;
  NJs(1) = 0;
  NJs(2) = 1;
  NJs(3) = -1;

  % compute derivatives of x and y wrt psi and eta
  xr = NJr*xcoords';
  yr = NJr*ycoords'; 
  xs = NJs*xcoords';
  ys = NJs*ycoords';

  % compute Jacobian
  Jinv = [ys, -yr; -xs, xr];
  jcob = xr*ys - xs*yr;
    
  % loop over 12 gauss points
  for g = 1:12
    % get parametric coordinates of current gauss point 'g'
    r = gr(g);
    s = gs(g);

    % shape functions
    N1 = r;
    N2 = s;
    N3 = 1-r-s;
    N = [ N1  0   N2  0   N3  0;
          0   N1  0   N2  0   N3];
        
    % get real coordinates of current gauss point
    gpcoords = N * elecoords;
    
    % approximate displacement at gauss point
    dis_approx = N * eledis;
    
    %% get "exact" displacement at gauss point
    % find the nearest node to the current gauss point
    % compute a distance vector
    dist_x = x_ref - gpcoords(1);
    dist_y = y_ref - gpcoords(2);
    dist_vec = sqrt(dist_x .^2 + dist_y .^2);
    
    % get nodeID with minimal distance
    [val nodeID] = min(dist_vec);
    
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
      p = [gpcoords(1) gpcoords(2)];

      % check, if current element 'aa' containts the gauss point 'p'
      if polygon_contains_point_2d (3,v,p);
        % 'aa' contains gauss point 'p'
        eleID_ref = aa;
      end;
    end;
    
    % get nodes of reference element
    elenodes_ref = node_ref(:,eleID_ref);
    
    % get nodal coordinates of reference element
    xcoords_ref = x_ref(elenodes_ref);
    ycoords_ref = y_ref(elenodes_ref);
    elecoords_ref = [ xcoords_ref(1) 
                      ycoords_ref(1)
                      xcoords_ref(2) 
                      ycoords_ref(2)
                      xcoords_ref(3) 
                      ycoords_ref(3)];
    
    % evaluate shape functions in reference mesh
    area_ref = det([[1 1 1]' xcoords_ref' ycoords_ref']) / 2;
    N_ref = zeros(2,6);
    % loop over nodes
    for b = 1:3     % Evaluate shape functions
      % load node coordinates
      xes_ref = xcoords_ref;
      yes_ref = ycoords_ref;

      % Get coordinates of area opposite node of concern
      xes_ref(b) = p(1); 
      yes_ref(b) = p(2);

      Larea_ref = det([[1 1 1]' xes_ref' yes_ref']) / 2;

      % Evaluate shape function for node 'b' (only in base DOFs)
      N_ref(1,2*b-1)  = N_ref(1,2*b-1)  + Larea_ref / area_ref;
      N_ref(2,2*b)    = N_ref(2,2*b)    + Larea_ref / area_ref;
    end;
    
    % get nodal displacements in reference mesh
    eledis_ref = [  dis_ref(2*elenodes_ref(1) - 1)
                    dis_ref(2*elenodes_ref(1))
                    dis_ref(2*elenodes_ref(2) - 1)
                    dis_ref(2*elenodes_ref(2))
                    dis_ref(2*elenodes_ref(3) - 1)
                    dis_ref(2*elenodes_ref(3))     ];
                  
    % compute "exact" displacement vector at gauss point
    dis_exact = N_ref * eledis_ref;
    % ------------------------------------------------------------------- %
    %% sum up the error
    % error at current gauss point
    dis_error = dis_approx - dis_exact
    
    % sqaure it
    dis_error2 = dis_error' * dis_error
    dis_exact2 = dis_exact' * dis_exact;
    
    % sum it up
    ele_error = ele_error + dis_error2 * wg(g) * jcob / 2;
    ele_solu = ele_solu + dis_exact2 * wg(g) * jcob / 2;
    % ------------------------------------------------------------------- %
  end;
  total_error = total_error + ele_error;
  total_solu = total_solu + ele_solu;
end;

% compute L2-norm
L2norm = sqrt(total_error / total_solu);
% ----------------------------------------------------------------------- %
%% print results
disp(['L2-norm of error in displacement field: ' num2str(L2norm)]);
% ----------------------------------------------------------------------- %