% DBCs

ubar = zeros(2,numnod);
dispbc = zeros(2,numnod);

for i=1:73
%     dispbc(1,i*19) = 1;
    dispbc(2,i*19) = 1;
end;

dispbc(1,19) = 1;
% dispbc(1,1) = 1;

