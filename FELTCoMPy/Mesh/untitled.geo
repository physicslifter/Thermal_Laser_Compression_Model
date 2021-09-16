// Gmsh project created on Wed Sep 15 12:24:22 2021
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, 6, 0, 1.0};
//+
Point(3) = {6, 6, 0, 1.0};
//+
Point(4) = {6, 4, 0, 1.0};
//+
Point(5) = {4, 4, 0, 1.0};
//+
Point(6) = {4, 2, 0, 1.0};
//+
Point(7) = {2, 2, 0, 1.0};
//+
Point(8) = {2, 0, 0, 1.0};
//+
Line(1) = {1, 8};
//+
Line(2) = {8, 7};
//+
Line(3) = {7, 6};
//+
Line(4) = {6, 5};
//+
Line(5) = {5, 4};
//+
Line(6) = {4, 3};
//+
Line(7) = {3, 2};
//+
Line(8) = {2, 1};
//+
Curve Loop(1) = {8, 1, 2, 3, 4, 5, 6, 7};
//+
Plane Surface(1) = {1};
//+
Physical Surface("plane_surface", 9) = {1};
//+
Physical Curve("initial_boundary", 10) = {8};
//+
Physical Curve("Fe-MGO_boundary1", 11) = {2};
//+
Physical Curve("Fe-MGO_boundary2", 12) = {4};
//+
Physical Curve("Fe-MGO_boundary3", 13) = {6};
//+
Physical Curve("top_wall_Fe", 14) = {7};
//+
Physical Curve("bottom_wall", 15) = {1, 3, 5};
