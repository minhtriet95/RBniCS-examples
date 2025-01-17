// Gmsh project created on Fri Feb  5 17:42:53 2021
//+
Point(1) = {0, 0, 1.5131, 1.0};
//+
Point(2) = {0, 0, -1.5131, 1.0};
//+
Point(3) = {0.2794, 0, -1.5131, 1.0};
//+
Point(4) = {0.2294, 0, -1.5131, 1.0};
//+
Point(5) = {0, 0.2794, -1.5131, 1.0};
//+
Point(6) = {0, 0.2294, -1.5131, 1.0};
//+
Point(7) = {0, 0.2294, 1.5131, 1.0};
//+
Point(8) = {0, 0.2794, 1.5131, 1.0};
//+
Point(9) = {0.2794, 0, 1.5131, 1.0};
//+
Point(10) = {0.2294, 0, 1.5131, 1.0};
//+
Circle(1) = {4, 2, 6};
//+
Circle(2) = {3, 2, 5};
//+
Circle(3) = {10, 1, 7};
//+
Circle(4) = {9, 1, 8};
//+
Line(5) = {5, 6};
//+
Line(6) = {4, 3};
//+
Line(7) = {8, 7};
//+
Line(8) = {10, 9};
//+
Curve Loop(1) = {1, -5, -2, -6};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {3, -7, -4, -8};
//+
Plane Surface(2) = {2};
//+
Line(9) = {4, 10};
//+
Line(10) = {7, 6};
//+
Line(11) = {5, 8};
//+
Line(12) = {9, 3};
//+
Curve Loop(3) = {1, -10, -3, -9};
//+
Surface(3) = {3};
//+
Curve Loop(4) = {2, 11, -4, 12};
//+
Surface(4) = {4};
//+
Point(11) = {0, 0, 1.7925, 1.0};
//+
Point(12) = {0, 0, 1.7425, 1.0};
//+
Point(13) = {0, 0, -1.7425, 1.0};
//+
Point(14) = {0, 0, -1.7925, 1.0};
//+
Circle(13) = {10, 1, 12};
//+
Line(15) = {12, 11};
//+
Line(16) = {13, 14};
//+
Circle(17) = {13, 2, 4};
//+
Circle(18) = {3, 2, 14};
//+
Circle(19) = {6, 2, 13};
//+
Circle(20) = {5, 2, 14};
//+
Circle(22) = {7, 1, 12};
//+
Curve Loop(5) = {14, -21, -4};
//+
Point(19) = {0, 0.05, 1.5131+(0.2794^2-0.05^2)^0.5, 1.0};
//+
Curve Loop(6) = {22, -13, 3};
//+
Surface(6) = {6};
//+
Curve Loop(7) = {22, 15, -21, 7};
//+
Curve Loop(8) = {14, -15, -13, 8};
//+
Curve Loop(9) = {19, 16, -20, 5};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {17, 6, 18, -16};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {20, -18, 2};
//+
Surface(11) = {11};
//+
Curve Loop(12) = {19, 17, 1};
//+
Surface(12) = {12};
//+
Curve Loop(13) = {10, -5, 11, 7};
//+
Plane Surface(13) = {13};
//+
Curve Loop(14) = {9, 8, 12, -6};
//+
Plane Surface(14) = {14};
//+
Surface Loop(1) = {13, 3, 14, 4, 1, 2};
//+
Surface Loop(2) = {1, 9, 12, 10, 11};
//+
Surface Loop(3) = {2, 8, 5, 7, 6};
//+
Point(15) = {0, 0, 1.8425, 1.0};
//+
Point(16) = {0, 0.05, 1.8425, 1.0};
//+
Point(17) = {0.05, 0, 1.8425, 1.0};
//+
Circle(23) = {16, 15, 17};
//+
Line(24) = {16, 15};
//+
Line(25) = {15, 17};
//+
Curve Loop(15) = {24, 25, -23};
//+
Plane Surface(15) = {15};
//+
Point(20) = {0.05, 0, 1.5131+(0.2794^2-0.05^2)^0.5, 1.0};
//+
Circle(26) = {19, 1, 11};
//+
Circle(14) = {20, 1, 11};
//+
Circle(27) = {9, 1, 20};
//+
Circle(28) = {8, 1, 19};
//+
Line(29) = {19, 16};
//+
Line(30) = {20, 17};
//+
Curve Loop(16) = {28, 26, -15, -22, -7};
//+
Plane Surface(16) = {16};
//+
Curve Loop(17) = {8, 27, 14, -15, -13};
//+
Plane Surface(17) = {17};
//+
Curve Loop(18) = {28, 26, -14, -27, 4};
//+
Line(31) = {11, 15};
//+
Curve Loop(19) = {29, 24, -31, -26};
//+
Plane Surface(18) = {19};
//+
Curve Loop(20) = {30, -25, -31, -14};
//+
Plane Surface(19) = {20};
//+
Circle(32) = {20, 11, 19};
//+
Curve Loop(21) = {27, 32, -28, -4};
//+
Surface(20) = {21};
//+
Curve Loop(22) = {30, -23, -29, -32};
//+
Surface(21) = {22};
//+
Curve Loop(23) = {32, 26, -14};
//+
Surface(22) = {23};
//+
Volume(1) = {1};
//+
Volume(2) = {2};
//+
Surface Loop(4) = {20, 17, 16, 6, 2, 22};
//+
Volume(3) = {4};
//+
Surface Loop(5) = {21, 19, 15, 18, 22};
//+
Volume(4) = {5};
//+
Physical Volume("Block 1", 1) = {1};
//+
Physical Volume("Block 2", 2) = {2};
//+
Physical Volume("Block 3", 3) = {3};
//+
Physical Volume("Block 4", 4) = {4};
//+
Physical Surface("Sideset 2", 2) = {9, 13, 16, 18};
//+
Physical Surface("Sideset 3", 3) = {10, 14, 17, 19};
//+
Physical Surface("Sideset 1", 1) = {15};
//+
Physical Surface("Sideset 4", 4) = {12, 3, 6};
