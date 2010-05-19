Point(1) = {0, 2, 0, 1e+022};
Point(2) = {0, -2, 0, 1e+022};
Point(3) = {16, -2, 0, 1e+022};
Point(4) = {16, 2, 0, 1e+022};
Line(1) = {1, 1};
Line(2) = {4, 3};
Line(3) = {3, 2};
Line(4) = {2, 1};
Line(5) = {1, 4};
Line Loop(7) = {2, 3, 4, 5};
Plane Surface(7) = {7};
Point(5) = {8.3, 2, 0};
Point(6) = {8.3, -2, 0};
Point(7) = {7.7, -2, 0};
Point(8) = {7.7, 2, 0};
Line(8) = {5, 6};
Line(9) = {7, 8};
Plane Surface(10) = {7};
Delete {
  Line{8, 9};
}
Delete {
  Point{5, 8, 7, 6};
}
