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
  %node = zeros(4,numele);
  numele = 0;
  for e=1:numsideele
      elemsideno = fscanf(fid, '%d',1);
      c1 = fscanf(fid, '%lf',1);
      if (c1 == 1) %a side
          c2 = fscanf(fid, '%lf',1);
          conn = fscanf(fid, '%lf',4);
      elseif (c1 == 2)
          numele = numele +1;
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