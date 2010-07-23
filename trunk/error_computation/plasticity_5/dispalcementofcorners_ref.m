% Compute the displacements of the bottom left and right corner.

% get mesh data
num_x = IFnldivx;
num_y = IFnldivy;

% bottom left corner
nodeID = num_x + 1;
dis_left = x_def(nodeID) - x_orig(nodeID)

% bottom right corner
nodeID = (num_y + 1) * (num_x + 1);
dis_right = x_def(nodeID) - x_orig(nodeID)