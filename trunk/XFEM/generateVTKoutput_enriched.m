disp('generate VTK data ...');

% call a subroutine to compute the displacement vectors for the nodes of the
% subelements
getVTKdata;

%% VTK-file for stresses and displacements
% create filename
filenameVTK = sprintf('VTK_files/VTK_%s_%d.vtk',filename_input_file,timestep);

% get an pointer to the output file
% options:
%   w         delete all contents and get writing acces
%   l         bit ordering scheme 'little endian'
fileID = fopen(filenameVTK,'w','l');

% write the file version and identifier
fileversionidentifier = sprintf('# vtk DataFile Version 2.0\n');
fprintf(fileID,fileversionidentifier);

% write the header
header = sprintf('Output data for %s\n', filename_input_file);
fprintf(fileID,header);

% define file format 'ASCII'
fileformat = sprintf('ASCII\n');
fprintf(fileID,fileformat);

% define nodal coordinates
dataset = sprintf('DATASET UNSTRUCTURED_GRID\nPOINTS %d  float\n',size(VTKcoords,1));
fprintf(fileID,dataset);
for i=1:size(VTKcoords,1)
  point = sprintf('%f %f 0\n',VTKcoords(i,1),VTKcoords(i,2));
  fprintf(fileID,point);
end;

% define element-node-connectivity
cells = sprintf('CELLS %d %d\n',size(VTKconn,1),4*size(VTKconn,1));
fprintf(fileID,cells);
for i=1:size(VTKconn,1)
  connectivity = sprintf('3 %d %d %d\n',VTKconn(i,1),VTKconn(i,2),VTKconn(i,3));
  fprintf(fileID,connectivity);
end;

% define cell types
celltypes = sprintf('CELL_TYPES %d\n',size(VTKconn,1));
fprintf(fileID,celltypes);
for i=1:size(VTKconn,1)
  fprintf(fileID,'5\n');
end;

% % define displacement field as a scalars with 3 components
% %   1   x-displacements
% %   2   y-displacements
% %   3   z-displacement (=0)
% fprintf(fileID,'POINT_DATA %d\n',size(VTKdis,1));
% fprintf(fileID,'SCALARS u(x) float 3\n');
% fprintf(fileID,'LOOKUP_TABLE default\n');
% for i=1:size(VTKdis,1);
%   fprintf(fileID,'%f %f 0.0\n',VTKdis(i,1),VTKdis(i,2));
% end;

fprintf(fileID,'POINT_DATA %d\n',size(VTKdis,1));
fprintf(fileID,'VECTORS u(x) float\n');
for i=1:size(VTKdis,1);
  fprintf(fileID,'%f %f 0.0\n',VTKdis(i,1),VTKdis(i,2));
end;

% define stresses in elements as CELL-data since stresses are constant in
% an element due to linear shape functions as SCALARS
%   1   xx-stress
%   2   yy-stress
%   3   xy-stress
%   von-Mises-stress
fprintf(fileID,'CELL_DATA %d\n',size(VTKstress,1));
fprintf(fileID,'SCALARS stress float 4\n');
fprintf(fileID,'LOOKUP_TABLE default\n');
for i=1:size(VTKstress,1);
  fprintf(fileID,'%f %f %f %f\n',VTKstress(i,1),VTKstress(i,2),VTKstress(i,3),VTKstress(i,4));
end;

% % define stresses in elements as CELL-data since stresses are constant in
% % an element due to linear shape functions as TENSORS
% fprintf(fileID,'CELL_DATA %d\n',size(VTKstress,1));
% fprintf(fileID,'TENSORS stress_tensor float\n');
% for i=1:size(VTKstress,1);
%   fprintf(fileID,'%f %f 0.0\n%f %f 0.0\n0.0 0.0 0.0\n\n',VTKstress(i,1),VTKstress(i,3),VTKstress(i,3),VTKstress(i,2));
% end;

% define skalar field, that shows, if a node is part of an element edge
% (=1) or not (=0). This is needed to show, that the cut elements are
% really separated into two parts, that slide along each other.
fprintf(fileID,'POINT_DATA %d\n',size(VTKcutgrid,1));
fprintf(fileID,'SCALARS cut_grid float 1\n');
fprintf(fileID,'LOOKUP_TABLE default\n');
for i=1:size(VTKcutgrid,1);
  fprintf(fileID,'%f \n',VTKcutgrid(i,1));
end;

% close file
status1 = fclose(fileID);
% ----------------------------------------------------------------------- %
%% VTK-file for stick-slip-zone
% create filename
filenameVTK = sprintf('VTK_files/VTK_%s_stickslip_%d.vtk',filename_input_file,timestep);

% get an pointer to the output file
% options:
%   w         delete all contents and get writing acces
%   l         bit ordering scheme 'little endian'
fileID = fopen(filenameVTK,'w','l');

% write the file version and identifier
fileversionidentifier = sprintf('# vtk DataFile Version 2.0\n');
fprintf(fileID,fileversionidentifier);

% write the header
header = sprintf('Evolution of stick-slip-zone for %s\n', filename_input_file);
fprintf(fileID,header);

% define file format 'ASCII'
fileformat = sprintf('ASCII\n');
fprintf(fileID,fileformat);

% define nodal coordinates
dataset = sprintf('DATASET UNSTRUCTURED_GRID\nPOINTS %d  float\n',size(VTKinterfacenodes,1));
fprintf(fileID,dataset);
for i=1:size(VTKinterfacenodes,1)
  point = sprintf('%f %f 0\n',VTKinterfacenodes(i,1),VTKinterfacenodes(i,2));
  fprintf(fileID,point);
end;

% define element-node-connectivity
cells = sprintf('CELLS %d %d\n',size(VTKinterfacestate,1),3*size(VTKinterfacestate,1));
fprintf(fileID,cells);
for i=1:2:2*size(VTKinterfacestate,1)
  connectivity = sprintf('2 %d %d \n',i-1,i);
  fprintf(fileID,connectivity);
end;

% define cell types
celltypes = sprintf('CELL_TYPES %d\n',size(VTKinterfacestate,1));
fprintf(fileID,celltypes);
for i=1:size(VTKinterfacestate,1)
  fprintf(fileID,'3\n');
end;

% define cell data
% represents stick or slip
% 0 ... stick
% 1 ... stick
fprintf(fileID,'CELL_DATA %d\n',size(VTKinterfacestate,1));
fprintf(fileID,'SCALARS stick_slip float 1\n');
fprintf(fileID,'LOOKUP_TABLE default\n');
for i=1:size(VTKinterfacestate,1);
  fprintf(fileID,'%f \n',VTKinterfacestate(i,1));
end;

% close file
status2 = fclose(fileID);
% ----------------------------------------------------------------------- %

if status1 == 0 && status2 == 0
  disp('VTK data generated successfully.');
end;