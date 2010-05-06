fid = fopen('mesh2.msh', 'r')
  c = fscanf(fid,'%s',1)
  numnod = fscanf(fid,'%d',1)
  x = zeros(1,numnod);
  y = zeros(1,numnod);
  
  for n=1:numnod
     nodeno = fscanf(fid,'%d',1);
     coord = fscanf(fid, '%lf',3);
     x(nodeno) = coord(1);
     y(nodeno) = coord(2);
  end
  
  c = fscanf(fid,'%s',1)
  c = fscanf(fid,'%s',1)
  
  numsideele = fscanf(fid,'%d',1)
   text = 'blabla3333'
  numele = 0;       % stores number of elements in the domain
  numboundele = 0;  % stores number of element edges on boundary of domain
  for e=1:numsideele
      elemsideno = fscanf(fid, '%d',1);
      c1 = fscanf(fid, '%lf',1);
      if (c1 == 1) %a side
          text = 'blabla22222'
          numboundele = numboundele + 1
          c2 = fscanf(fid, '%lf',3);
          conn = fscanf(fid, '%lf',2);
          boundary_nodes(1:2,numboundele);
      elseif (c1 == 2)
          text = 'blabla'
          numele = numele + 1;
          c2 = fscanf(fid, '%lf',1); 
          c3 = fscanf(fid, '%lf',1); 
          c4 = fscanf(fid, '%lf',1);
          conn = fscanf(fid, '%lf',3);
          node(1:3,numele) = conn';    
      elseif (c1 == 3)
          numele = numele +1;
          c2 = fscanf(fid, '%lf',1); 
          c3 = fscanf(fid, '%lf',1); 
          c4 = fscanf(fid, '%lf',1);
          conn = fscanf(fid, '%lf',4);
          node(1:4,numele) = conn';
      end
      
  end
  
fclose(fid)