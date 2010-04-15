function [xint_new,edge_ids] = extend_segment(elemi, p1, p2, cuteleminfo);

%routine to extend the line of a segment to determine other intersection
%point
global X Y CONN;
tol = 0.001;

if (cuteleminfo.nb_int > 1)
    error('in extend_segment1','the nb of intersections is greater than one');
else
    if (cuteleminfo.edge_ids(1) == 1) %check sides 2 and 3
        %side 2
        q1 = [X(CONN(2,elemi)) Y(CONN(2,elemi))]; 
        q2 = [X(CONN(3,elemi)) Y(CONN(3,elemi))];
        [ flag1, p ] = lines_exp_int_2d ( p1, p2, q1, q2 );
        if (flag1 == 1) %check to make sure it's actually on the segment
            [ u, v ] = segment_contains_point_2d ( q1, q2, p );
            if ( u >= -0.000000001 & u <= 1.000000001 & v < tol)
              xint_new(1,1:2) = cuteleminfo.xint;
              xint_new(2,1:2) = p;
              edge_ids = [1 2];
              return
            end
        end
        %side 3
        q1 = [X(CONN(3,elemi)) Y(CONN(3,elemi))]; 
        q2 = [X(CONN(1,elemi)) Y(CONN(1,elemi))];
        [ flag2, p ] = lines_exp_int_2d ( p1, p2, q1, q2 );
        if (flag2 == 1) 
            [ u, v ] = segment_contains_point_2d ( q1, q2, p );
            if ( u >= -0.000000001 & u <= 1.000000001 & v < tol)
               xint_new(1,1:2) = cuteleminfo.xint;
               xint_new(2,1:2) = p;
               edge_ids = [1 3];
               return
            end
        end
        
        %if you get to this point, no other intersection has been detected
        flag1
        flag2
        xint_new(1,1:2) = cuteleminfo.xint;
        xint_new(2,1:2) = cuteleminfo.xint;
        edge_ids = [1 1]
        %error('in extend segment2','No intersection detected');
        
    elseif (cuteleminfo.edge_ids(1) == 2) %check sides 1 and 3
        %side 1
        q1 = [X(CONN(1,elemi)) Y(CONN(1,elemi))]; 
        q2 = [X(CONN(2,elemi)) Y(CONN(2,elemi))];
        [ flag1, p ] = lines_exp_int_2d ( p1, p2, q1, q2 );
        if (flag1 == 1) %check to make sure it's actually on the segment
            [ u, v ] = segment_contains_point_2d ( q1, q2, p );
            if ( u >= -0.000000001 & u <= 1.000000001 & v < tol)
              xint_new(1,1:2) = p;
              xint_new(2,1:2) = cuteleminfo.xint;
              edge_ids = [1 2];
              return
            end
        end
        %side 3
        q1 = [X(CONN(3,elemi)) Y(CONN(3,elemi))]; 
        q2 = [X(CONN(1,elemi)) Y(CONN(1,elemi))];
        [ flag2, p ] = lines_exp_int_2d ( p1, p2, q1, q2 );
        if (flag2 == 1) 
            [ u, v ] = segment_contains_point_2d ( q1, q2, p );
            if ( u >= -0.000000001 & u <= 1.00000001 & v < tol)
               xint_new(1,1:2) = cuteleminfo.xint;
               xint_new(2,1:2) = p;
               edge_ids = [2 3];
               return
            end
        end
        
        %if you get to this point, no other intersection has been detected
        %given the fact that at least one point was intersected, this means
        %a degenerate case - likely an intersection with a node, so there's
        %no need to extend.  simply return the same point which will kick
        %out later
        flag1
        flag2
        xint_new(1,1:2) = cuteleminfo.xint;
        xint_new(2,1:2) = cuteleminfo.xint;
        edge_ids = [2 2]
        %error('in extend segment3','No intersection detected');
        
    elseif (cuteleminfo.edge_ids(1) == 3) %check sides 1 and 2
        %side 1
        q1 = [X(CONN(1,elemi)) Y(CONN(1,elemi))]; 
        q2 = [X(CONN(2,elemi)) Y(CONN(2,elemi))];
        [ flag1, p ] = lines_exp_int_2d ( p1, p2, q1, q2 );
        if (flag1 == 1) %check to make sure it's actually on the segment
            [ u, v ] = segment_contains_point_2d ( q1, q2, p );
            if ( u >= -0.000000001 & u <= 1.000000001 & v < tol)
              xint_new(1,1:2) = p;
              xint_new(2,1:2) = cuteleminfo.xint;
              edge_ids = [1 3];
              return
            end
        end
        
        %side 2
        q1 = [X(CONN(2,elemi)) Y(CONN(2,elemi))]; 
        q2 = [X(CONN(3,elemi)) Y(CONN(3,elemi))];
        [ flag2, p ] = lines_exp_int_2d ( p1, p2, q1, q2 );
        if (flag2 == 1) %check to make sure it's actually on the segment
            [ u, v ] = segment_contains_point_2d ( q1, q2, p );
            if ( u >= -0.000000001 & u <= 1.000000001 & v < tol)
              xint_new(1,1:2) = p;
              xint_new(2,1:2) = cuteleminfo.xint;
              edge_ids = [2 3];
              return
            end
        end
        
        %if you get to this point, no other intersection has been detected
        flag1
        flag2
        xint_new(1,1:2) = cuteleminfo.xint;
        xint_new(2,1:2) = cuteleminfo.xint;
        edge_ids = [3 3]
        %error('in extend segment4','No intersection detected');
        
    else
        cuteleminfo.edge_ids
        error('in extend_segment5','the edge_ids do not make sense');
    end
end