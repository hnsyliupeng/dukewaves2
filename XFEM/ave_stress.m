% Program to calculate the L2 norm of the tractions across the boundary 
% in the "cstruct" convergence test.  Calculation of - <sigma> dot n

% Only one boundary exists, at x = 5 y = variable

% Initial L2norm

total_errorx = 0;
total_errory = 0;
total_tracx = 0;
total_tracy = 0;

for i = 1:numele  % loop over elements
    
    if (cutlist(i) ~= 0)     % If the element is cut by an interface
        
        % Number of sub-interfaces
        
        num_sub_int = size(INT_INTERFACE(i).pairings,1);
        
        for j = 1:num_sub_int   % Loop over sub interface
            i;
            
            % Get normal to interface
            
            %[norm] = get_norm(i,j);
            norm = [1 0];

            young = 1000;
            pr = 0.0;
            fac = young/(1 - (pr)^2);
            D = fac*[1.0, pr, 0;
                    pr, 1.0, 0.0;
                    0, 0, (1.-pr)/2 ];

      
            % get coordinates of element nodes 
            for m=1:3
                je = node(m,i); xe(m) = x(je); ye(m) = y(je);
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
                           
            % element displacement vector
            for m=1:3
                % base degrees of freedom
                m1 = id_eqns(node(m,i),1);
                m2 = id_eqns(node(m,i),2); 
                dispj(m*2-1) = fdisp(m1);
                dispj(m*2) = fdisp(m2);
            end
            
            %base stress
            stress1 = (D*BJ*dispj')';

            % element displacement vector
            for m=1:3
                % extra degrees of freedom
                m1 = id_eqns(node(m,i),3);
                m2 = id_eqns(node(m,i),4); 
                dispj(m*2-1) = fdisp(m1);
                dispj(m*2) = fdisp(m2);
            end

            %enriched stress
            stress2 = (D*BJ*dispj')';
            
            % average stress
            sigma = stress1 + stress2/2;
            
      
            % Get y coordinates of sub interface
            
            ycoord = [INT_INTERFACE(i).shared(2*j-1,2)
                        INT_INTERFACE(i).shared(2*j,2)];
            
            j;
            
            sig = [sigma(1) sigma(3);
                   sigma(3) sigma(2)]
                 
            % Dot sigma and norm
            
            trac_temp = - sig*norm';
            
            % Grab the x tractions
            
            tracx = [trac_temp(1) trac_temp(1)];  % Set it up to plot well
            tracy = [trac_temp(2) trac_temp(2)];     
            
            % plot
            
            figure(1)
            plot(ycoord,tracx)
            hold on
            figure(2)
            plot(ycoord,tracy)
            hold on
            
            % Now, compute the L2norm of the error over the sub interface
            
            % Two point gauss quadrature
            
            % jacobian of segment to global
            he = abs(ycoord(1)-ycoord(2));
            seg_jcob = he/2;
            
            gauss = [-sqrt(3)/3 sqrt(3)/3];
            weights = [1 1];
            
            segment_errx = 0;
            segment_erry = 0;
            traction_intx = 0;
            traction_inty = 0;

            for g = 1:2
                      
                % Get real coordinate at Gauss point
                yn = 0.5*(1-gauss(g))*ycoord(1)+0.5*(1+gauss(g))*ycoord(2);
                
                % Get value of function at that point               
                anal(1) = -2.0625*yn;
                anal(2) = 3/32*(4-yn^2);
                
                % calc = analytical - error!
                errx =  trac_temp(1) - anal(1);
                erry =  trac_temp(2) - anal(2);
                
                errx2 = errx^2;
                erry2 = erry^2;
                tracx2 = trac_temp(1)^2;
                tracy2 = trac_temp(2)^2;
                
                
                % assemble error squared over segment

                segment_errx = segment_errx + errx2*weights(g)*seg_jcob;
                segment_erry = segment_erry + erry2*weights(g)*seg_jcob;
                traction_intx = traction_intx + tracx2*weights(g)*seg_jcob;
                traction_inty = traction_inty + tracy2*weights(g)*seg_jcob;
                
            end
            
            total_errorx = total_errorx + segment_errx;
            total_errory = total_errory + segment_erry;
            total_tracx = total_tracx + traction_intx; 
            total_tracy = total_tracy + traction_inty;  
        end
    end
end

L2trac = sqrt(total_errorx + total_errory)/sqrt(total_tracx + total_tracy);

% Plot analytical solns

coord = linspace(-2,2);
real_x = (-1)*(16 - 5)*coord/(16/3);
real_y = 3/32*(4-coord.^2);

figure(1)
plot(coord,real_x,'--')
title('Analytical and Numerical Traction Fields')
xlabel('Distance Along interface')
ylabel('x-component of traction')

figure(2)
plot(coord,real_y,'--')
title('Analytical and Numerical Traction Fields')
xlabel('Distance Along interface')
ylabel('y-component of traction')
