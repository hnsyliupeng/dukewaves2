function ExportFieldGmsh(x,y,z,node,u)

fid = fopen('Field.dat','w');
%header information
  fprintf(fid, '$MeshFormat \n');
  fprintf(fid, '2.0 0 8 \n');
  fprintf(fid, '$EndMeshFormat \n');
  
%node data  
  numnod = size(x,1);
  fprintf(fid, '$Nodes \n');
  fprintf(fid, '%d \n', numnod);
  
  for i=1:numnod
    fprintf(fid, '%d %7.4f %7.4f %7.4f \n', i, x(i), y(i), z(i));
  end
  fprintf(fid, '$EndNodes \n');
  
%Element data
  numele = size(node,2);
  fprintf(fid, '$Elements \n');
  fprintf(fid, '%d \n', numele);
  
  for e=1:numele
     fprintf(fid, '%d 4 2 99 2 ', e); %4 here indicates the element is a tet
     for j=1:size(node(:,e),1)
       fprintf(fid, '%d ', node(j,e));
     end
     fprintf(fid, ' \n');
  end
  fprintf(fid, '$EndElements \n');
  
%Field Information
  fprintf(fid, '$NodeData \n');
  fprintf(fid, '1 \n');  
  fprintf(fid, '\"Approximate solution u^h\" \n');
  fprintf(fid, '1 \n');
  fprintf(fid, '0.0 \n'); %the time value
  fprintf(fid, '3 \n'); %three integer tags
  fprintf(fid, '0 \n'); %the time step
  
  fprintf(fid, '1 \n'); %1 component (for a scalar field)
  fprintf(fid, '%d \n', size(u,1) );  %give the number of values
  for j=1:size(u,1)
    fprintf(fid, '%d %8.4f \n', j, u(j));
  end
  
  fprintf(fid, '$EndNodeData \n');
  
fclose(fid);
