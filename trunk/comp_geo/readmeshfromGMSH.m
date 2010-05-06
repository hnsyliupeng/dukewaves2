% readmeshfromGMSH.m
%
% reads mesh from GMSH-format '*.msh' into variables, that can be used in
% xfem-computation. It works only for triangular elements.
%
% Nodes, that are use for drawing the geometry, but which are not part of
% the mesh, are called 'non physical' in the following. They have to be
% rejected. Else, the stiffness matrix is singular.
%

% build filename
filename_msh_file = [IFfilename_msh_file '.msh'];
filename_msh_file = fullfile(pwd,'GMSHmeshes',filename_msh_file);

fid = fopen(filename_msh_file, 'r');
c = fscanf(fid,'%s',6);
numnod = fscanf(fid,'%d',1);
x = zeros(1,numnod);
y = zeros(1,numnod);

for n=1:numnod
    numnod = fscanf(fid,'%d',1);
    coord = fscanf(fid, '%lf',3);
    x(numnod) = coord(1);
    y(numnod) = coord(2);
end

c = fscanf(fid,'%s',1);
c = fscanf(fid,'%s',1);


nonphysicalnodes = 1:numnod;

numsideele = fscanf(fid,'%d',1);
%node = zeros(3,numele);
numele = 0;
numboundele = 0;
%   node = [];
%   boundary_nodes = [];
for e=1:numsideele
    elemsideno = fscanf(fid, '%d',1);
    c1 = fscanf(fid, '%lf',1);
    if (c1 == 1) %a side
        c2 = fscanf(fid, '%lf',4);
        conn = fscanf(fid, '%lf',2);
        if conn(1) ~= conn(2)     % repress equal nodes
            numboundele = numboundele + 1;
            boundary_nodes(1:2,numboundele) = conn;
        end;
    elseif (c1 == 2)
        numele = numele +1;
        c2 = fscanf(fid, '%lf',4); 
%           c3 = fscanf(fid, '%lf',1); 
%           c4 = fscanf(fid, '%lf',1);
%           c5 = fscanf(fid, '%lf',1);
        conn = fscanf(fid, '%lf',3);
        node(1:3,numele) = conn';    
        for i=1:3
            nonphysicalnodes(conn(i))=0;
        end;
    elseif (c1 == 3)
%           numele = numele +1;
        c2 = fscanf(fid, '%lf',8); 
%           c3 = fscanf(fid, '%lf',1); 
%           c4 = fscanf(fid, '%lf',1);
%           conn = fscanf(fid, '%lf',4);
%           node(1:4,numele) = conn';
    elseif (c1 == 15)
        c2 = fscanf(fid, '%lf', 5);
    else
        error('MATLAB:readmesh:UnvalidMeshType',...
            'Mesh consists elements, that are no triangles.');
    end

end
  
fclose(fid);

index = 0;
nonphysnodevec = [0];
for nodeid = 1:numnod
    if nonphysicalnodes(nodeid) ~= 0
%         index = index + 1;
%         x_temp(index) = x(nodeid);
%         y_temp(index) = y(nodeid);
        nonphysnodevec = [nonphysnodevec nonphysicalnodes(nodeid)];    
    end;
end;

nonphysnodevec = nonphysnodevec(2:end);

% for j=1:3
%     for k=1:numele
%         for i=nonphysnodevec            
%             if node(j,k) > i
%                 node(j,k) = node(j,k) - 1;
%             end;
%         end;
%     end;
% end;

% x = x_temp;
% y = y_temp;
% node = node_temp;
% numnod = index;
% numele = index2;

% print some data
text1 = ['numnod:       ' num2str(numnod)];
text2 = ['numele:       ' num2str(numele)];
text3 = ['numboundele:  ' num2str(numboundele)];
disp(text1);
disp(text2);
disp(text3);

% clear all variables, that are not necessary to describe a mesh
clear c c1 c2 conn coord e elemsideno fid n numsideele x_temp y_temp...
    nodeid nonphysicalnodes index i j k node_temp text1 text2; 
