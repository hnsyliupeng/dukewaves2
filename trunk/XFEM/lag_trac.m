% Program to calculate the L2 norm of the tractions across the boundary in
% the cstruct convergence test from the penalty contribution.  Calculation
% of \beta \int{[[u]]}

% Initialize L2norm

total_errorx = 0;
total_errory = 0;
total_tracx = 0;
total_tracy = 0;


cnt = 0;
for i = 1:numele  % loop over elements
    
    if (cutlist(i) ~= 0)     % If the element is cut by an interface
        
        cnt = cnt + 1;
        
        % Number of sub-interfaces
        
        num_sub_int = size(INT_INTERFACE(i).pairings,1);
        
        for j = 1:num_sub_int   % Loop over sub interface

            % end points of intersection - direction doesn't matter - this is for the
            % segment jacobian calculation

            p1 = [INT_INTERFACE(i).shared(2*j-1,1) INT_INTERFACE(i).shared(2*j-1,2)];
            p2 = [INT_INTERFACE(i).shared(2*j,1) INT_INTERFACE(i).shared(2*j,2)];
            
            % Get y coordinates of sub interface
            
            ycoord = [p1(2),p2(2)];
            
            % Grab the x tractions
            
            tracx = [-fdisp(old_size + 2*cnt - 1) -fdisp(old_size + cnt*2 - 1)];  
            tracy = [-fdisp(old_size + 2*cnt) -fdisp(old_size + 2*cnt)];     
            
            % plot
           
            figure(1)
            plot(ycoord,tracx)
            hold on
            figure(2)
            plot(ycoord,tracy)
            hold on
            
            % jacobian of segment to global
            he = sqrt((p1(1)-p2(1))^2 + (p1(2)-p2(2))^2);
            seg_jcob = he/2;

            % Gauss points on segments
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
               anal(1) = -1*(16-4.4)/(16/3)*yn;
               anal(2) = 3/32*(4-yn^2);
               
               % calc = analytical - error!
               errx =  tracx(1) - anal(1);
               erry =  tracy(2) - anal(2);
               
               errx2 = errx^2;
               erry2 = erry^2;
               tracx2 = tracx(1)^2;
               tracy2 = tracy(2)^2;
                
                
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
real_x = (-1)*(16 - 4.4)*coord/(16/3);
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

