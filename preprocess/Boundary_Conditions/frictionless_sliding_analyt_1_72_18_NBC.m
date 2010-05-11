%NBCs

% maximum load in x-direction
maxload = -1/72;

force = zeros(2,numnod);

for i=0:72
    force(2,19*i+1)=i*maxload;
end;

% force(2,1369)=-1;
