% DBCs

% fixed nodes
% fixed in x-direction
for i=1:19
    dispbc(1,i) = 1;
%     dispbc(2,i) = 1;
end;

% fixed in y-direction
dispbc(2,9) = 1;