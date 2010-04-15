function [nb_sub_elts, sub_elt_ids] = slice_2d_element_by_segment(elemi, cuteleminfo);

global X Y CONN;
%first compare the element number to that in the cuteleminfo
%if they're not the same then that indicates the latter needs to be
%generated

%tolerance to check against a segment intersecting a node of the element
tol = 1.0e-06;

debug_slice2d = 0;

maxelemid = size(CONN,2);
maxnodeid = size(X,2);

if (elemi == cuteleminfo.elemno)
      xint = cuteleminfo.xint;
      %the above should always have three points, below four (note this is
      %completely arbitrary - I'm not using a positive or negative
      %unless there is a degenerate case - one of the intersection points
      %overlaps a node
      %'side' of the interface
      if (cuteleminfo.edge_ids == [1 2])
          %first check overlap of xint with element nodes 1 and 3
          check1 = check_overlap(xint(1,1), xint(1,2), X(CONN(1,elemi)), Y(CONN(1,elemi)),tol);
          check3 = check_overlap(xint(2,1), xint(2,2), X(CONN(3,elemi)), Y(CONN(3,elemi)),tol);
          if (check1)
              %only insert the second intersection point into X and Y
              X(1,maxnodeid+1) = xint(2,1); Y(1,maxnodeid+1) = xint(2,2);
              xp_abo = [X(CONN(2,elemi)) xint(2,1) X(CONN(1,elemi))];
              yp_abo = [Y(CONN(2,elemi)) xint(2,2) Y(CONN(1,elemi))];
              xptoX_abo = [CONN(2,elemi) maxnodeid+1 CONN(1,elemi)];
              xp_bel = [X(CONN(1,elemi)) X(CONN(3,elemi)) xint(2,1)];
              yp_bel = [Y(CONN(1,elemi)) Y(CONN(3,elemi)) xint(2,2)];
              xptoX_bel = [CONN(1,elemi) CONN(3,elemi) maxnodeid+1];
          elseif (check3)
              %only insert the first intersection point into X and Y
              X(1,maxnodeid+1) = xint(1,1); Y(1,maxnodeid+1) = xint(1,2);
              xp_abo = [X(CONN(2,elemi)) X(CONN(3,elemi)) xint(1,1) ];
              yp_abo = [Y(CONN(2,elemi)) Y(CONN(3,elemi)) xint(1,2) ];
              xptoX_abo = [CONN(2,elemi) CONN(3,elemi) maxnodeid+1];
              xp_bel = [X(CONN(1,elemi)) X(CONN(3,elemi)) xint(1,1) ];
              yp_bel = [Y(CONN(1,elemi)) Y(CONN(3,elemi)) xint(1,2) ];
              xptoX_bel = [CONN(1,elemi) CONN(3,elemi) maxnodeid+1 ];
          else
             %insert the points into the list of nodes
             X(1,maxnodeid+1) = xint(1,1); X(1,maxnodeid+2) = xint(2,1);
             Y(1,maxnodeid+1) = xint(1,2); Y(1,maxnodeid+2) = xint(2,2);
             xp_abo = [X(CONN(2,elemi)) xint(1,1) xint(2,1)];
             yp_abo = [Y(CONN(2,elemi)) xint(1,2) xint(2,2)];
             xptoX_abo = [CONN(2,elemi) maxnodeid+1 maxnodeid+2];
             xp_bel = [X(CONN(1,elemi)) X(CONN(3,elemi)) xint(1,1) xint(2,1)];
             yp_bel = [Y(CONN(1,elemi)) Y(CONN(3,elemi)) xint(1,2) xint(2,2)];
             xptoX_bel = [CONN(1,elemi) CONN(3,elemi) maxnodeid+1 maxnodeid+2];
          end
      elseif (cuteleminfo.edge_ids == [1 3])
          %first check overlap of xint with element nodes 2 and 3
          check2 = check_overlap(xint(1,1), xint(1,2), X(CONN(2,elemi)), Y(CONN(2,elemi)),tol);
          check3 = check_overlap(xint(2,1), xint(2,2), X(CONN(3,elemi)), Y(CONN(3,elemi)),tol);
          if (check2)
             %only insert the second point
             X(1,maxnodeid+1) = xint(2,1); Y(1,maxnodeid+1) = xint(2,2); 
             xp_abo = [X(CONN(1,elemi)) X(CONN(2,elemi)) xint(2,1)];
             yp_abo = [Y(CONN(1,elemi)) Y(CONN(2,elemi)) xint(2,2)];
             xptoX_abo = [CONN(1,elemi) CONN(2,elemi) maxnodeid+1 ];
             xp_bel = [X(CONN(2,elemi)) X(CONN(3,elemi)) xint(2,1)];
             yp_bel = [Y(CONN(2,elemi)) Y(CONN(3,elemi)) xint(2,2)];
             xptoX_bel = [CONN(2,elemi) CONN(3,elemi) maxnodeid+1];
          elseif (check3)
             %only insert the first point
             X(1,maxnodeid+1) = xint(1,1); Y(1,maxnodeid+1) = xint(1,2);
             xp_abo = [X(CONN(1,elemi)) X(CONN(3,elemi)) xint(1,1) ];
             yp_abo = [Y(CONN(1,elemi)) Y(CONN(3,elemi)) xint(1,2) ];
             xptoX_abo = [CONN(1,elemi) CONN(3,elemi) maxnodeid+1];
             xp_bel = [X(CONN(2,elemi)) X(CONN(3,elemi)) xint(1,1) ];
             yp_bel = [Y(CONN(2,elemi)) Y(CONN(3,elemi)) xint(1,2) ];
             xptoX_bel = [CONN(2,elemi) CONN(3,elemi) maxnodeid+1 ];
          else
             %insert the points into the list of nodes
             X(1,maxnodeid+1) = xint(1,1); X(1,maxnodeid+2) = xint(2,1);
             Y(1,maxnodeid+1) = xint(1,2); Y(1,maxnodeid+2) = xint(2,2);
             xp_abo = [X(CONN(1,elemi)) xint(1,1) xint(2,1)];
             yp_abo = [Y(CONN(1,elemi)) xint(1,2) xint(2,2)];
             xptoX_abo = [CONN(1,elemi) maxnodeid+1 maxnodeid+2];
             xp_bel = [X(CONN(2,elemi)) X(CONN(3,elemi)) xint(1,1) xint(2,1)];
             yp_bel = [Y(CONN(2,elemi)) Y(CONN(3,elemi)) xint(1,2) xint(2,2)];
             xptoX_bel = [CONN(2,elemi) CONN(3,elemi) maxnodeid+1 maxnodeid+2];
          end
      elseif (cuteleminfo.edge_ids == [2 3])
          %first check overlap of xint with element nodes 1 and 2
          check1 = check_overlap(xint(2,1), xint(2,2), X(CONN(1,elemi)), Y(CONN(1,elemi)),tol);
          check2 = check_overlap(xint(1,1), xint(1,2), X(CONN(2,elemi)), Y(CONN(2,elemi)),tol);
          if (check1)
             %only insert the first point
             X(1,maxnodeid+1) = xint(1,1); Y(1,maxnodeid+1) = xint(1,2);
             xp_abo = [X(CONN(3,elemi)) X(CONN(1,elemi)) xint(1,1) ];
             yp_abo = [Y(CONN(3,elemi)) Y(CONN(1,elemi)) xint(1,2) ];
             xptoX_abo = [CONN(3,elemi) CONN(1,elemi) maxnodeid+1 ];
             xp_bel = [X(CONN(1,elemi)) X(CONN(2,elemi)) xint(1,1) ];
             yp_bel = [Y(CONN(1,elemi)) Y(CONN(2,elemi)) xint(1,2) ];
             xptoX_bel = [CONN(1,elemi) CONN(2,elemi) maxnodeid+1 ];
          elseif (check2)
             %only insert the second point
             X(1,maxnodeid+1) = xint(2,1); Y(1,maxnodeid+1) = xint(2,2);
             xp_abo = [X(CONN(3,elemi)) xint(2,1) X(CONN(2,elemi))];
             yp_abo = [Y(CONN(3,elemi)) xint(2,2) Y(CONN(2,elemi))];
             xptoX_abo = [CONN(3,elemi) maxnodeid+1 CONN(2,elemi)];
             xp_bel = [X(CONN(1,elemi)) X(CONN(2,elemi)) xint(2,1)];
             yp_bel = [Y(CONN(1,elemi)) Y(CONN(2,elemi)) xint(2,2)];
             xptoX_bel = [CONN(1,elemi) CONN(2,elemi) maxnodeid+1 ];
          else
             %insert the points into the list of nodes
             X(1,maxnodeid+1) = xint(1,1); X(1,maxnodeid+2) = xint(2,1);
             Y(1,maxnodeid+1) = xint(1,2); Y(1,maxnodeid+2) = xint(2,2);
             xp_abo = [X(CONN(3,elemi)) xint(1,1) xint(2,1)];
             yp_abo = [Y(CONN(3,elemi)) xint(1,2) xint(2,2)];
             xptoX_abo = [CONN(3,elemi) maxnodeid+1 maxnodeid+2];
             xp_bel = [X(CONN(1,elemi)) X(CONN(2,elemi)) xint(1,1) xint(2,1)];
             yp_bel = [Y(CONN(1,elemi)) Y(CONN(2,elemi)) xint(1,2) xint(2,2)];
             xptoX_bel = [CONN(1,elemi) CONN(2,elemi) maxnodeid+1 maxnodeid+2];
          end
      else
          elemi
          cuteleminfo.edge_ids
          error('slice_2d_element_by_segment1','The edge_ids are not correct');
      end
      %generate the new subelements
      tri_abo = delaunay(xp_abo,yp_abo);
      tri_bel = delaunay(xp_bel, yp_bel);
      %add them to the connectivity, above then below
      if (size(tri_abo,1) == 1)
          CONN(1:3,maxelemid+1) = xptoX_abo(tri_abo(1,1:3));
      else
          tri_abo
          error('slice_2d_element_by_segment2','The above triangulation failed');
      end
      if (size(tri_bel,1) == 1)
          CONN(1:3,maxelemid+2) = xptoX_bel(tri_bel(1,1:3));
      elseif (size(tri_bel,1) == 2)
          CONN(1:3,maxelemid+2) = xptoX_bel(tri_bel(1,1:3));
          CONN(1:3,maxelemid+3) = xptoX_bel(tri_bel(2,1:3));
      else
          tri_bel
          error('slice_2d_element_by_segment3','The below triangulation failed');
      end
      nb_sub_elts = size(tri_abo,1)+size(tri_bel,1);
      if (nb_sub_elts == 2)
         sub_elt_ids = [maxelemid+1 maxelemid+2];
      elseif (nb_sub_elts == 3)
         sub_elt_ids = [maxelemid+1 maxelemid+2 maxelemid+3];
      else
          error('slice_2d_element_by_segment3b','The nb_sub_elts is not 2 or 3');
      end
      %plot the elements if debug
      if (debug_slice2d)
         figure(1)
         hold on
         for j=1:size(sub_elt_ids,2);
            e = sub_elt_ids(j);
            plot(X(CONN(1:3,e)),Y(CONN(1:3,e)),'g')
            plot([X(CONN(3,e)) X(CONN(1,e))],[Y(CONN(3,e)) Y(CONN(1,e))],'g')
         end
         drawnow
      end
else %generate it - this is for a subelement 
    p1 = cuteleminfo.xint(1,1:2);
    p2 = cuteleminfo.xint(2,1:2);
    nb_int = 0;
    for j=1:2 %loop over edges to determine intersection
       q1 = [X(CONN(j,elemi)) Y(CONN(j,elemi))];
       q2 = [X(CONN(j+1,elemi)) Y(CONN(j+1,elemi))];
       [flag, pint] = segments_int_2d(p1,p2,q1,q2);
       if (flag)
          nb_int = nb_int+1;
          xint(nb_int,1:2) = pint;
          edge_ids(nb_int) = j;
       end
    end
    q1 = [X(CONN(3,elemi)) Y(CONN(3,elemi))]; 
    q2 = [X(CONN(1,elemi)) Y(CONN(1,elemi))];
    [flag, pint] = segments_int_2d(p1,p2,q1,q2);
    if (flag)
        nb_int = nb_int+1;
        xint(nb_int,1:2) = pint;
        edge_ids(nb_int) = 3;
    end
    
    if (nb_int == 3) %has to have been a degenerate case
      nb_sub_elts = 1;
      sub_elt_ids = -1;
      return 
    end
    
    if (nb_int ==0 )
        nb_sub_elts = 1;
        sub_elt_ids = -1;
        return
    elseif (nb_int == 1)
        %extend it first
        cuteleminfo_sub = struct('elemno',elemi,'nb_int',nb_int,'xint',xint,'edge_ids',edge_ids);        
        [xint_new,edge_ids_new] = extend_segment(elemi, p1, p2, cuteleminfo_sub);
        clear xint edge_ids;
        xint = xint_new; edge_ids = edge_ids_new;
        
        %check if the two intersection points aren't identical
        check0 = check_overlap(xint(1,1), xint(1,2), xint(2,1), xint(2,2), tol);
        if (check0)
          nb_sub_elts = 1;
          sub_elt_ids = -1;
         return
        end
        
        if (edge_ids == [1 2])
          %first check overlap of xint with element nodes 1 and 3
          check1 = check_overlap(xint(1,1), xint(1,2), X(CONN(1,elemi)), Y(CONN(1,elemi)),tol);
          check3 = check_overlap(xint(2,1), xint(2,2), X(CONN(3,elemi)), Y(CONN(3,elemi)),tol);
          if (check1)
              %only insert the second intersection point into X and Y
              X(1,maxnodeid+1) = xint(2,1); Y(1,maxnodeid+1) = xint(2,2);
              xp_abo = [X(CONN(2,elemi)) xint(2,1) X(CONN(1,elemi))];
              yp_abo = [Y(CONN(2,elemi)) xint(2,2) Y(CONN(1,elemi))];
              xptoX_abo = [CONN(2,elemi) maxnodeid+1 CONN(1,elemi)];
              xp_bel = [X(CONN(1,elemi)) X(CONN(3,elemi)) xint(2,1)];
              yp_bel = [Y(CONN(1,elemi)) Y(CONN(3,elemi)) xint(2,2)];
              xptoX_bel = [CONN(1,elemi) CONN(3,elemi) maxnodeid+1];
          elseif (check3)
              %only insert the first intersection point into X and Y
              X(1,maxnodeid+1) = xint(1,1); Y(1,maxnodeid+1) = xint(1,2);
              xp_abo = [X(CONN(2,elemi)) X(CONN(3,elemi)) xint(1,1) ];
              yp_abo = [Y(CONN(2,elemi)) Y(CONN(3,elemi)) xint(1,2) ];
              xptoX_abo = [CONN(2,elemi) CONN(3,elemi) maxnodeid+1];
              xp_bel = [X(CONN(1,elemi)) X(CONN(3,elemi)) xint(1,1) ];
              yp_bel = [Y(CONN(1,elemi)) Y(CONN(3,elemi)) xint(1,2) ];
              xptoX_bel = [CONN(1,elemi) CONN(3,elemi) maxnodeid+1 ];
          else
             %insert the points into the list of nodes
             X(1,maxnodeid+1) = xint(1,1); X(1,maxnodeid+2) = xint(2,1);
             Y(1,maxnodeid+1) = xint(1,2); Y(1,maxnodeid+2) = xint(2,2);
             xp_abo = [X(CONN(2,elemi)) xint(1,1) xint(2,1)];
             yp_abo = [Y(CONN(2,elemi)) xint(1,2) xint(2,2)];
             xptoX_abo = [CONN(2,elemi) maxnodeid+1 maxnodeid+2];
             xp_bel = [X(CONN(1,elemi)) X(CONN(3,elemi)) xint(1,1) xint(2,1)];
             yp_bel = [Y(CONN(1,elemi)) Y(CONN(3,elemi)) xint(1,2) xint(2,2)];
             xptoX_bel = [CONN(1,elemi) CONN(3,elemi) maxnodeid+1 maxnodeid+2];
          end
      elseif (edge_ids == [1 3])
          %first check overlap of xint with element nodes 2 and 3
          check2 = check_overlap(xint(1,1), xint(1,2), X(CONN(2,elemi)), Y(CONN(2,elemi)),tol);
          check3 = check_overlap(xint(2,1), xint(2,2), X(CONN(3,elemi)), Y(CONN(3,elemi)),tol);
          if (check2)
             %only insert the second point
             X(1,maxnodeid+1) = xint(2,1); Y(1,maxnodeid+1) = xint(2,2); 
             xp_abo = [X(CONN(1,elemi)) X(CONN(2,elemi)) xint(2,1)];
             yp_abo = [Y(CONN(1,elemi)) Y(CONN(2,elemi)) xint(2,2)];
             xptoX_abo = [CONN(1,elemi) CONN(2,elemi) maxnodeid+1 ];
             xp_bel = [X(CONN(2,elemi)) X(CONN(3,elemi)) xint(2,1)];
             yp_bel = [Y(CONN(2,elemi)) Y(CONN(3,elemi)) xint(2,2)];
             xptoX_bel = [CONN(2,elemi) CONN(3,elemi) maxnodeid+1];
          elseif (check3)
             %only insert the first point
             X(1,maxnodeid+1) = xint(1,1); Y(1,maxnodeid+1) = xint(1,2);
             xp_abo = [X(CONN(1,elemi)) X(CONN(3,elemi)) xint(1,1) ];
             yp_abo = [Y(CONN(1,elemi)) Y(CONN(3,elemi)) xint(1,2) ];
             xptoX_abo = [CONN(1,elemi) CONN(3,elemi) maxnodeid+1];
             xp_bel = [X(CONN(2,elemi)) X(CONN(3,elemi)) xint(1,1) ];
             yp_bel = [Y(CONN(2,elemi)) Y(CONN(3,elemi)) xint(1,2) ];
             xptoX_bel = [CONN(2,elemi) CONN(3,elemi) maxnodeid+1 ];
          else
             %insert the points into the list of nodes
             X(1,maxnodeid+1) = xint(1,1); X(1,maxnodeid+2) = xint(2,1);
             Y(1,maxnodeid+1) = xint(1,2); Y(1,maxnodeid+2) = xint(2,2);
             xp_abo = [X(CONN(1,elemi)) xint(1,1) xint(2,1)];
             yp_abo = [Y(CONN(1,elemi)) xint(1,2) xint(2,2)];
             xptoX_abo = [CONN(1,elemi) maxnodeid+1 maxnodeid+2];
             xp_bel = [X(CONN(2,elemi)) X(CONN(3,elemi)) xint(1,1) xint(2,1)];
             yp_bel = [Y(CONN(2,elemi)) Y(CONN(3,elemi)) xint(1,2) xint(2,2)];
             xptoX_bel = [CONN(2,elemi) CONN(3,elemi) maxnodeid+1 maxnodeid+2];
          end
      elseif (edge_ids == [2 3])
          %first check overlap of xint with element nodes 1 and 2
          check1 = check_overlap(xint(2,1), xint(2,2), X(CONN(1,elemi)), Y(CONN(1,elemi)),tol);
          check2 = check_overlap(xint(1,1), xint(1,2), X(CONN(2,elemi)), Y(CONN(2,elemi)),tol);
          if (check1)
             %only insert the first point
             X(1,maxnodeid+1) = xint(1,1); Y(1,maxnodeid+1) = xint(1,2);
             xp_abo = [X(CONN(3,elemi)) X(CONN(1,elemi)) xint(1,1) ];
             yp_abo = [Y(CONN(3,elemi)) Y(CONN(1,elemi)) xint(1,2) ];
             xptoX_abo = [CONN(3,elemi) CONN(1,elemi) maxnodeid+1 ];
             xp_bel = [X(CONN(1,elemi)) X(CONN(2,elemi)) xint(1,1) ];
             yp_bel = [Y(CONN(1,elemi)) Y(CONN(2,elemi)) xint(1,2) ];
             xptoX_bel = [CONN(1,elemi) CONN(2,elemi) maxnodeid+1 ];
          elseif (check2)
             %only insert the second point
             X(1,maxnodeid+1) = xint(2,1); Y(1,maxnodeid+1) = xint(2,2);
             xp_abo = [X(CONN(3,elemi)) xint(2,1) X(CONN(2,elemi))];
             yp_abo = [Y(CONN(3,elemi)) xint(2,2) Y(CONN(2,elemi))];
             xptoX_abo = [CONN(3,elemi) maxnodeid+1 CONN(2,elemi)];
             xp_bel = [X(CONN(1,elemi)) X(CONN(2,elemi)) xint(2,1)];
             yp_bel = [Y(CONN(1,elemi)) Y(CONN(2,elemi)) xint(2,2)];
             xptoX_bel = [CONN(1,elemi) CONN(2,elemi) maxnodeid+1 ];
          else
             %insert the points into the list of nodes
             X(1,maxnodeid+1) = xint(1,1); X(1,maxnodeid+2) = xint(2,1);
             Y(1,maxnodeid+1) = xint(1,2); Y(1,maxnodeid+2) = xint(2,2);
             xp_abo = [X(CONN(3,elemi)) xint(1,1) xint(2,1)];
             yp_abo = [Y(CONN(3,elemi)) xint(1,2) xint(2,2)];
             xptoX_abo = [CONN(3,elemi) maxnodeid+1 maxnodeid+2];
             xp_bel = [X(CONN(1,elemi)) X(CONN(2,elemi)) xint(1,1) xint(2,1)];
             yp_bel = [Y(CONN(1,elemi)) Y(CONN(2,elemi)) xint(1,2) xint(2,2)];
             xptoX_bel = [CONN(1,elemi) CONN(2,elemi) maxnodeid+1 maxnodeid+2];
          end
      else
          elemi
          cuteleminfo.edge_ids
          error('slice_2d_element_by_segment1','The edge_ids are not correct');
      end
      %generate the new subelements
      tri_abo = delaunay(xp_abo,yp_abo);
      tri_bel = delaunay(xp_bel, yp_bel);
      %add them to the connectivity, above then below
      if (size(tri_abo,1) == 1)
          CONN(1:3,maxelemid+1) = xptoX_abo(tri_abo(1,1:3));
      else
          tri_abo
          error('slice_2d_element_by_segment2','The above triangulation failed');
      end
      if (size(tri_bel,1) == 1)
          CONN(1:3,maxelemid+2) = xptoX_bel(tri_bel(1,1:3));
      elseif (size(tri_bel,1) == 2)
          CONN(1:3,maxelemid+2) = xptoX_bel(tri_bel(1,1:3));
          CONN(1:3,maxelemid+3) = xptoX_bel(tri_bel(2,1:3));
      else
          tri_bel
          error('slice_2d_element_by_segment3','The below triangulation failed');
      end
      nb_sub_elts = size(tri_abo,1)+size(tri_bel,1);
      if (nb_sub_elts == 2)
         sub_elt_ids = [maxelemid+1 maxelemid+2];
      elseif (nb_sub_elts == 3)
         sub_elt_ids = [maxelemid+1 maxelemid+2 maxelemid+3];
      else
          error('slice_2d_element_by_segment3b','The nb_sub_elts is not 2 or 3');
      end
      %plot the elements if debug
      if (debug_slice2d)
         figure(1)
         hold on
         for j=1:size(sub_elt_ids,2);
            e = sub_elt_ids(j);
            plot(X(CONN(1:3,e)),Y(CONN(1:3,e)),'g')
            plot([X(CONN(3,e)) X(CONN(1,e))],[Y(CONN(3,e)) Y(CONN(1,e))],'g')
         end
         drawnow
      end
      return
    elseif (nb_int == 2)
        
       %nb_int
       %xint
       
       %check if the two intersection points aren't identical
       check0 = check_overlap(xint(1,1), xint(1,2), xint(2,1), xint(2,2), tol);
       if (check0)
         nb_sub_elts = 1;
         sub_elt_ids = -1;
        return
       end
        
       if (edge_ids == [1 2])
          %first check overlap of xint with element nodes 1 and 3
          check1 = check_overlap(xint(1,1), xint(1,2), X(CONN(1,elemi)), Y(CONN(1,elemi)),tol);
          check3 = check_overlap(xint(2,1), xint(2,2), X(CONN(3,elemi)), Y(CONN(3,elemi)),tol);
          if (check1)
              %only insert the second intersection point into X and Y
              X(1,maxnodeid+1) = xint(2,1); Y(1,maxnodeid+1) = xint(2,2);
              xp_abo = [X(CONN(2,elemi)) xint(2,1) X(CONN(1,elemi))];
              yp_abo = [Y(CONN(2,elemi)) xint(2,2) Y(CONN(1,elemi))];
              xptoX_abo = [CONN(2,elemi) maxnodeid+1 CONN(1,elemi)];
              xp_bel = [X(CONN(1,elemi)) X(CONN(3,elemi)) xint(2,1)];
              yp_bel = [Y(CONN(1,elemi)) Y(CONN(3,elemi)) xint(2,2)];
              xptoX_bel = [CONN(1,elemi) CONN(3,elemi) maxnodeid+1];
          elseif (check3)
              %only insert the first intersection point into X and Y
              X(1,maxnodeid+1) = xint(1,1); Y(1,maxnodeid+1) = xint(1,2);
              xp_abo = [X(CONN(2,elemi)) X(CONN(3,elemi)) xint(1,1) ];
              yp_abo = [Y(CONN(2,elemi)) Y(CONN(3,elemi)) xint(1,2) ];
              xptoX_abo = [CONN(2,elemi) CONN(3,elemi) maxnodeid+1];
              xp_bel = [X(CONN(1,elemi)) X(CONN(3,elemi)) xint(1,1) ];
              yp_bel = [Y(CONN(1,elemi)) Y(CONN(3,elemi)) xint(1,2) ];
              xptoX_bel = [CONN(1,elemi) CONN(3,elemi) maxnodeid+1 ];
          else
             %insert the points into the list of nodes
             X(1,maxnodeid+1) = xint(1,1); X(1,maxnodeid+2) = xint(2,1);
             Y(1,maxnodeid+1) = xint(1,2); Y(1,maxnodeid+2) = xint(2,2);
             xp_abo = [X(CONN(2,elemi)) xint(1,1) xint(2,1)];
             yp_abo = [Y(CONN(2,elemi)) xint(1,2) xint(2,2)];
             xptoX_abo = [CONN(2,elemi) maxnodeid+1 maxnodeid+2];
             xp_bel = [X(CONN(1,elemi)) X(CONN(3,elemi)) xint(1,1) xint(2,1)];
             yp_bel = [Y(CONN(1,elemi)) Y(CONN(3,elemi)) xint(1,2) xint(2,2)];
             xptoX_bel = [CONN(1,elemi) CONN(3,elemi) maxnodeid+1 maxnodeid+2];
          end
      elseif (edge_ids == [1 3])
          %first check overlap of xint with element nodes 2 and 3
          check2 = check_overlap(xint(1,1), xint(1,2), X(CONN(2,elemi)), Y(CONN(2,elemi)),tol);
          check3 = check_overlap(xint(2,1), xint(2,2), X(CONN(3,elemi)), Y(CONN(3,elemi)),tol);
          if (check2)
             %only insert the second point
             X(1,maxnodeid+1) = xint(2,1); Y(1,maxnodeid+1) = xint(2,2); 
             xp_abo = [X(CONN(1,elemi)) X(CONN(2,elemi)) xint(2,1)];
             yp_abo = [Y(CONN(1,elemi)) Y(CONN(2,elemi)) xint(2,2)];
             xptoX_abo = [CONN(1,elemi) CONN(2,elemi) maxnodeid+1 ];
             xp_bel = [X(CONN(2,elemi)) X(CONN(3,elemi)) xint(2,1)];
             yp_bel = [Y(CONN(2,elemi)) Y(CONN(3,elemi)) xint(2,2)];
             xptoX_bel = [CONN(2,elemi) CONN(3,elemi) maxnodeid+1];
          elseif (check3)
             %only insert the first point
             X(1,maxnodeid+1) = xint(1,1); Y(1,maxnodeid+1) = xint(1,2);
             xp_abo = [X(CONN(1,elemi)) X(CONN(3,elemi)) xint(1,1) ];
             yp_abo = [Y(CONN(1,elemi)) Y(CONN(3,elemi)) xint(1,2) ];
             xptoX_abo = [CONN(1,elemi) CONN(3,elemi) maxnodeid+1];
             xp_bel = [X(CONN(2,elemi)) X(CONN(3,elemi)) xint(1,1) ];
             yp_bel = [Y(CONN(2,elemi)) Y(CONN(3,elemi)) xint(1,2) ];
             xptoX_bel = [CONN(2,elemi) CONN(3,elemi) maxnodeid+1 ];
          else
             %insert the points into the list of nodes
             X(1,maxnodeid+1) = xint(1,1); X(1,maxnodeid+2) = xint(2,1);
             Y(1,maxnodeid+1) = xint(1,2); Y(1,maxnodeid+2) = xint(2,2);
             xp_abo = [X(CONN(1,elemi)) xint(1,1) xint(2,1)];
             yp_abo = [Y(CONN(1,elemi)) xint(1,2) xint(2,2)];
             xptoX_abo = [CONN(1,elemi) maxnodeid+1 maxnodeid+2];
             xp_bel = [X(CONN(2,elemi)) X(CONN(3,elemi)) xint(1,1) xint(2,1)];
             yp_bel = [Y(CONN(2,elemi)) Y(CONN(3,elemi)) xint(1,2) xint(2,2)];
             xptoX_bel = [CONN(2,elemi) CONN(3,elemi) maxnodeid+1 maxnodeid+2];
          end
      elseif (edge_ids == [2 3])
          %first check overlap of xint with element nodes 1 and 2
          check1 = check_overlap(xint(2,1), xint(2,2), X(CONN(1,elemi)), Y(CONN(1,elemi)),tol);
          check2 = check_overlap(xint(1,1), xint(1,2), X(CONN(2,elemi)), Y(CONN(2,elemi)),tol);
          if (check1)
             %only insert the first point
             X(1,maxnodeid+1) = xint(1,1); Y(1,maxnodeid+1) = xint(1,2);
             xp_abo = [X(CONN(3,elemi)) X(CONN(1,elemi)) xint(1,1) ];
             yp_abo = [Y(CONN(3,elemi)) Y(CONN(1,elemi)) xint(1,2) ];
             xptoX_abo = [CONN(3,elemi) CONN(1,elemi) maxnodeid+1 ];
             xp_bel = [X(CONN(1,elemi)) X(CONN(2,elemi)) xint(1,1) ];
             yp_bel = [Y(CONN(1,elemi)) Y(CONN(2,elemi)) xint(1,2) ];
             xptoX_bel = [CONN(1,elemi) CONN(2,elemi) maxnodeid+1 ];
          elseif (check2)
             %only insert the second point
             X(1,maxnodeid+1) = xint(2,1); Y(1,maxnodeid+1) = xint(2,2);
             xp_abo = [X(CONN(3,elemi)) xint(2,1) X(CONN(2,elemi))];
             yp_abo = [Y(CONN(3,elemi)) xint(2,2) Y(CONN(2,elemi))];
             xptoX_abo = [CONN(3,elemi) maxnodeid+1 CONN(2,elemi)];
             xp_bel = [X(CONN(1,elemi)) X(CONN(2,elemi)) xint(2,1)];
             yp_bel = [Y(CONN(1,elemi)) Y(CONN(2,elemi)) xint(2,2)];
             xptoX_bel = [CONN(1,elemi) CONN(2,elemi) maxnodeid+1 ];
          else
             %insert the points into the list of nodes
             X(1,maxnodeid+1) = xint(1,1); X(1,maxnodeid+2) = xint(2,1);
             Y(1,maxnodeid+1) = xint(1,2); Y(1,maxnodeid+2) = xint(2,2);
             xp_abo = [X(CONN(3,elemi)) xint(1,1) xint(2,1)];
             yp_abo = [Y(CONN(3,elemi)) xint(1,2) xint(2,2)];
             xptoX_abo = [CONN(3,elemi) maxnodeid+1 maxnodeid+2];
             xp_bel = [X(CONN(1,elemi)) X(CONN(2,elemi)) xint(1,1) xint(2,1)];
             yp_bel = [Y(CONN(1,elemi)) Y(CONN(2,elemi)) xint(1,2) xint(2,2)];
             xptoX_bel = [CONN(1,elemi) CONN(2,elemi) maxnodeid+1 maxnodeid+2];
          end
      else
          elemi
          cuteleminfo.edge_ids
          error('slice_2d_element_by_segment1','The edge_ids are not correct');
      end
      %generate the new subelements
      tri_abo = delaunay(xp_abo,yp_abo);
      %xp_bel
      %yp_bel
      tri_bel = delaunay(xp_bel, yp_bel);
      %add them to the connectivity, above then below
      if (size(tri_abo,1) == 1)
          CONN(1:3,maxelemid+1) = xptoX_abo(tri_abo(1,1:3));
      else
          tri_abo
          error('slice_2d_element_by_segment2','The above triangulation failed');
      end
      if (size(tri_bel,1) == 1)
          CONN(1:3,maxelemid+2) = xptoX_bel(tri_bel(1,1:3));
      elseif (size(tri_bel,1) == 2)
          CONN(1:3,maxelemid+2) = xptoX_bel(tri_bel(1,1:3));
          CONN(1:3,maxelemid+3) = xptoX_bel(tri_bel(2,1:3));
      else
          tri_bel
          error('slice_2d_element_by_segment3','The below triangulation failed');
      end
      nb_sub_elts = size(tri_abo,1)+size(tri_bel,1);
      if (nb_sub_elts == 2)
         sub_elt_ids = [maxelemid+1 maxelemid+2];
      elseif (nb_sub_elts == 3)
         sub_elt_ids = [maxelemid+1 maxelemid+2 maxelemid+3];
      else
          error('slice_2d_element_by_segment3b','The nb_sub_elts is not 2 or 3');
      end
      %plot the elements if debug
      if (debug_slice2d)
         figure(1)
         hold on
         for j=1:size(sub_elt_ids,2);
            e = sub_elt_ids(j);
            plot(X(CONN(1:3,e)),Y(CONN(1:3,e)),'g')
            plot([X(CONN(3,e)) X(CONN(1,e))],[Y(CONN(3,e)) Y(CONN(1,e))],'g')
         end
         drawnow
      end
      return
    %else
    %    elemi
    %    error('slice_2d_element_by_segment7','Only one intersection found')
    end
end


