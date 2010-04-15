%partition the elements based on intersections with grain boundaries
%first, get some general information

%associate with each segment in the vor_d a list of cutelements
%for each element cut by a segment, store the number of intersections and
%locations

debug_part2d = 1;

%loop over segments
for seg_id = 1:size(vx,2)
  nb_cut_elems = 0;
  p1 = [vx(1,seg_id) vy(1,seg_id)]; p2 = [vx(2,seg_id) vy(2,seg_id)];
  for e=1:numele
     if (cutlist(e) > 0) 
         nb_int = 0; %the number of intersections between the current segment/element pair
         clear xint edge_ids;
         for j=1:2 %loop over edges to determine intersection
           q1 = [x(node(j,e)) y(node(j,e))]; q2 = [x(node(j+1,e)) y(node(j+1,e))];
           [flag, pint] = segments_int_2d(p1,p2,q1,q2);
             if (flag)
               nb_int = nb_int+1;
               xint(nb_int,1:2) = pint;
               edge_ids(nb_int) = j;
             end
         end
         q1 = [x(node(3,e)) y(node(3,e))]; q2 = [x(node(1,e)) y(node(1,e))];
         [flag, pint] = segments_int_2d(p1,p2,q1,q2);
         if (flag)
           nb_int = nb_int+1;
           xint(nb_int,1:2) = pint;
           edge_ids(nb_int) = 3;
         end
         if (nb_int > 0)
            nb_cut_elems = nb_cut_elems+1;
            cuteleminfo = struct('elemno',e,'nb_int',nb_int,'xint',xint,'edge_ids',edge_ids);
            seg_cut_info(seg_id,nb_cut_elems) = cuteleminfo;
         end        
     end %cut check
  end %loop over cut elements
  %there may be some segments that don't intersect the mesh at all
  if ( nb_cut_elems == 0) %initialize to test against
      cuteleminfo = struct('elemno',-1,'nb_int',0,'xint',[],'edge_ids',[]);
      seg_cut_info(seg_id,1) = cuteleminfo;
  end
end %loop

clear cuteleminfo;

if (debug_part2d)
   figure(3)
   plot(vx(1:2,1),vy(1:2,1),'r')
   hold on
   for i=1:size(seg_cut_info(1,:),2)
       e = seg_cut_info(1,i).elemno;
       if (e ~= -1)
          plot(x(node(1:3,e)),y(node(1:3,e)))
          plot([x(node(3,e)) x(node(1,e))],[y(node(3,e)) y(node(1,e))])
       end
       drawnow
   end
end

%stop

%initialize partition info
for e=1:numele
    eleminfo = struct('nb_subelts',1,'subelemids',-1); %use -1 to indicate the original element
    eleminfo_arr(e) = eleminfo;
end
global ELEMINFO_ARR;
ELEMINFO_ARR = eleminfo_arr;

%now slice recursively, setting up partition info along the way
for seg_id = 1:size(vx,2)
    seg_id
%for seg_id = 1:5
      for i=1:size(seg_cut_info(seg_id,:),2) %loop over elements cut by the segment
          elemi = seg_cut_info(seg_id,i).elemno;
          if (elemi>0) %tests against segments that are initialized to -1 (no intersections)
             if (seg_cut_info(seg_id,i).nb_int == 2)
                 cuteleminfo = seg_cut_info(seg_id,i);
                 recursive_slice_2d(elemi, cuteleminfo)
             elseif (seg_cut_info(seg_id,i).nb_int == 1)
                 %extend the segment into the element as necessary
                 %if (ELEMINFO_ARR(elemi).nb_subelts == 1)
                    cuteleminfo = seg_cut_info(seg_id,i);
                    p1 = [vx(1,seg_id) vy(1,seg_id)]; p2 = [vx(2,seg_id) vy(2,seg_id)];
                    [xint_new,edge_ids] = extend_segment(elemi, p1, p2, cuteleminfo);
                    %xint_new
                    %edge_ids
                    cuteleminfo.xint = xint_new;
                    cuteleminfo.edge_ids = edge_ids;
                    cuteleminfo.nb_int = 2;
                    %call the recursive 2d partition routine
                    recursive_slice_2d(elemi, cuteleminfo)
                 %else
                    % recursive_slice_2d(elemi, cuteleminfo)
                 %end
             else
                 seg_cut_info(seg_id,i).nb_int
                 error('partition2d','the number of cuts is not 1 or 2');
             end
          end
      end
end

