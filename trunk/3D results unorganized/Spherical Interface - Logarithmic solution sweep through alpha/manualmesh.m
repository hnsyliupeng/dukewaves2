function [x,y,z,node,numtet,numnod] = manualmesh(zdiv)
numdivz = zdiv;
numcube = numdivz;
numtet = 6*numcube;
numnod = 4*(numdivz+1);
node = zeros(4,numtet);
x = zeros(numnod,1);
y = zeros(numnod,1);
z = zeros(numnod,1);


nodefirstcube =[3,7,6,5,7,1; 
                6,6,1,7,4,8; 
                2,3,7,8,3,4; 
                1,1,5,1,1,7];

for j=1:(numdivz +1)
    if(j==1)
        z(1:4*j) = (j-1)*(1/numdivz);
        x(1:4*j) = [0;1;1;0];
        y(1:4*j) = [0;0;1;1];
    else
        z((4*(j-1)+1):4*j) = ((j-1)/numdivz);
        x((4*(j-1)+1):4*j) = [0;1;1;0];
        y((4*(j-1)+1):4*j) = [0;0;1;1];
    end
end

for ncube=1:numcube
    if(ncube==1)
        for i=1:4
            for j=1:6
                node(i,j) = node(i,j) + nodefirstcube(i,j);
            end
        end
    else
        for i=1:4
            for j=(6*(ncube-1)+1):6*ncube
                node(i,j) = node(i,j) + nodefirstcube(i,(j-6*(ncube-1))) + 4*(ncube-1);
            end
        end
    end
end
