% NBC

% mesh data
num_x = 21;%11;%21;%41;%61;%81;
num_y = 27;%13;%27;%51;%77;%101;

p = 0.001;         % load per length
loadlength = 4;   % length, that is loaded
load = p * loadlength / num_x; % nodal load

tractionnodes = num_x * (num_y + 1) + 1: num_x * (num_y + 1) + num_x;

for nodeID = tractionnodes
  force(1,nodeID) = load;
end;
force(1,tractionnodes(1)) = force(1,tractionnodes(1)) / 2;
force(1,tractionnodes(end)) = force(1,tractionnodes(end)) / 2;

clear tractionnodes;