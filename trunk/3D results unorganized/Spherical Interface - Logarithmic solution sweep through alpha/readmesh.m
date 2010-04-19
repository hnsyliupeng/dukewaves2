fprintf('Entering readmesh routine \n')
fid = fopen('Charlength-0.25.msh', 'r');
%fid = fopen('mesh1graded.msh', 'r')
  c = fscanf(fid,'%s',1);
  c = fscanf(fid,'%d',3);
  c = fscanf(fid,'%s',1);
  c = fscanf(fid,'%s',1);
  numnod = fscanf(fid,'%d',1);
  x = zeros(numnod,1);
  y = zeros(numnod,1);
  z = zeros(numnod,1);
  
  fprintf('Number of nodes = %d \n', numnod)
 
  for n=1:numnod
     nodeno = fscanf(fid,'%d',1);
     coord = fscanf(fid, '%lf',3);
     x(nodeno) = coord(1);
     y(nodeno) = coord(2);
     z(nodeno) = coord(3);
  end
  
  c = fscanf(fid,'%s',1);
  c = fscanf(fid,'%s',1);
  
  numsideele = fscanf(fid,'%d',1);
  %node = zeros(4,numele);
  numele = 0;
  numtet = 0;
  for e=1:numsideele
      elemsideno = fscanf(fid, '%d',1);
      c1 = fscanf(fid, '%lf',1);
      if (c1 == 15)
          c2 = fscanf(fid, '%lf',1);
          c3 = fscanf(fid, '%lf',1);
          c4 = fscanf(fid, '%lf',1);
          c5 = fscanf(fid, '%lf',1);
          conn = fscanf(fid, '%lf',1);
      elseif (c1 == 1) %a side
          c2 = fscanf(fid, '%lf',1);
          c3 = fscanf(fid, '%lf',1);
          c4 = fscanf(fid, '%lf',1);
          c5 = fscanf(fid, '%lf',1);
          conn = fscanf(fid, '%lf',2);
      elseif (c1 == 2)
          numele = numele +1;
          c2 = fscanf(fid, '%lf',1); 
          c3 = fscanf(fid, '%lf',1); 
          c4 = fscanf(fid, '%lf',1);
          c5 = fscanf(fid, '%lf',1);
          conn = fscanf(fid, '%lf',3);
%           node(1:3,numele) = conn'    
      elseif (c1 == 3)
          numele = numele +1;
          c2 = fscanf(fid, '%lf',1); 
          c3 = fscanf(fid, '%lf',1); 
          c4 = fscanf(fid, '%lf',1);
          c5 = fscanf(fid, '%lf',1);
          conn = fscanf(fid, '%lf',4);
%           node(1:4,numele) = conn';
      elseif (c1 == 4)
          numtet = numtet +1;
          numele = numele +1;
          c2 = fscanf(fid, '%lf',1); 
          c3 = fscanf(fid, '%lf',1); 
          c4 = fscanf(fid, '%lf',1);
          c5 = fscanf(fid, '%lf',1);
          conn = fscanf(fid, '%lf',4);
          node(1:4,numtet) = conn';
          
      end
  end
  fprintf('Read mesh, exiting \n');
fclose(fid);