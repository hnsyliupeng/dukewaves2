Point(1) = {0, 1, 0, 1e+022};
Point(2) = {5, 1, 0, 1e+022};
Point(3) = {5, -1, 0, 1e+022};
Point(4) = {0, -1, 0, 1e+022};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(6) = {1, 2, 3, 4};
Plane Surface(6) = {6};