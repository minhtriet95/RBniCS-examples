//+
Point(1) = {-1, -1, 0, 1.0};
//+
Point(2) = {-1, 1, 0, 1.0};
//+
Point(3) = {1, 1, 0, 1.0};
//+
Point(4) = {1, -1, 0, 1.0};
//+
Point(5) = {0, 0, 0, 0};
//+
Point(6) = {0, 0.5, 0, 1.0};
//+
Point(7) = {0, -0.5, 0, 1.0};
//+
Line(1) = {1, 4};
//+
Line(2) = {4, 3};
//+
Line(3) = {3, 2};
//+
Line(4) = {2, 1};
//+
Point(8) = {-0.5, 0, 0, 1.0};
//+
Point(9) = {0.5, 0, 0, 1.0};
//+
Line(9) = {2, 8};
//+
Line(10) = {8, 1};
//+
Line(11) = {1, 7};
//+
Line(12) = {7, 4};
//+
Line(13) = {4, 9};
//+
Line(14) = {9, 3};
//+
Line(15) = {3, 6};
//+
Line(16) = {6, 2};
//+
Point(10) = {0.6, 0.6, 0, 1.0};
//+
Point(11) = {0.6, -0.6, 0, 1.0};
//+
Point(12) = {-0.6, -0.6, 0, 1.0};
//+
Point(13) = {-0.6, 0.6, 0, 1.0};
//+
Point(14) = {0.125^0.5, 0.125^0.5, 0, 1.0};
//+
Point(15) = {0.125^0.5, -0.125^0.5, 0, 1.0};
//+
Point(16) = {-0.125^0.5, -0.125^0.5, 0, 1.0};
//+
Point(17) = {-0.125^0.5, 0.125^0.5, 0, 1.0};
//+
Circle(17) = {7, 5, 15};
//+
Circle(18) = {15, 5, 9};
//+
Circle(19) = {9, 5, 14};
//+
Circle(20) = {14, 5, 6};
//+
Circle(21) = {6, 5, 17};
//+
Circle(22) = {17, 5, 8};
//+
Circle(23) = {8, 5, 16};
//+
Circle(24) = {16, 5, 7};
//+
Line(25) = {13, 8};
//+
Line(26) = {13, 6};
//+
Line(27) = {6, 10};
//+
Line(28) = {10, 9};
//+
Line(29) = {9, 11};
//+
Line(30) = {11, 7};
//+
Line(31) = {7, 12};
//+
Line(32) = {12, 8};
//+
Line(33) = {16, 12};
//+
Line(34) = {12, 1};
//+
Line(35) = {17, 13};
//+
Line(36) = {13, 2};
//+
Line(37) = {14, 10};
//+
Line(38) = {10, 3};
//+
Line(39) = {15, 11};
//+
Line(40) = {11, 4};
//+
Curve Loop(1) = {22, 23, 24, 17, 18, 19, 20, 21};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {24, 31, -33};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {33, 32, 23};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {25, -22, 35};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {21, 35, 26};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {27, -37, 20};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {19, 37, 28};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {18, 29, -39};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {17, 39, 30};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {11, 31, 34};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {34, -10, -32};
//+
Plane Surface(11) = {11};
//+
Curve Loop(12) = {9, -25, 36};
//+
Plane Surface(12) = {12};
//+
Curve Loop(13) = {26, 16, -36};
//+
Plane Surface(13) = {13};
//+
Curve Loop(14) = {15, 27, 38};
//+
Plane Surface(14) = {14};
//+
Curve Loop(15) = {28, 14, -38};
//+
Plane Surface(15) = {15};
//+
Curve Loop(16) = {13, 29, 40};
//+
Plane Surface(16) = {16};
//+
Curve Loop(17) = {12, -40, 30};
//+
Plane Surface(17) = {17};
//+
Curve Loop(18) = {11, 12, -1};
//+
Plane Surface(18) = {18};
//+
Curve Loop(19) = {10, -4, 9};
//+
Plane Surface(19) = {19};
//+
Curve Loop(20) = {16, -3, 15};
//+
Plane Surface(20) = {20};
//+
Curve Loop(21) = {14, -2, 13};
//+
Plane Surface(21) = {21};
//+
Physical Surface("Domain 1", 1) = {1};
//+
Physical Surface("Domain 2", 2) = {2};
//+
Physical Surface("Domain 3", 3) = {3};
//+
Physical Surface("Domain 4", 4) = {4};
//+
Physical Surface("Domain 5", 5) = {5};
//+
Physical Surface("Domain 6", 6) = {6};
//+
Physical Surface("Domain 7", 7) = {7};
//+
Physical Surface("Domain 8", 8) = {8};
//+
Physical Surface("Domain 9", 9) = {9};
//+
Physical Surface("Domain 10", 10) = {10};
//+
Physical Surface("Domain 11", 11) = {11};
//+
Physical Surface("Domain 12", 12) = {12};
//+
Physical Surface("Domain 13", 13) = {13};
//+
Physical Surface("Domain 14", 14) = {14};
//+
Physical Surface("Domain 15", 15) = {15};
//+
Physical Surface("Domain 16", 16) = {16};
//+
Physical Surface("Domain 17", 17) = {17};
//+
Physical Surface("Domain 18", 18) = {18};
//+
Physical Surface("Domain 19", 19) = {19};
//+
Physical Surface("Domain 20", 20) = {20};
//+
Physical Surface("Domain 21", 21) = {21};
//+
Physical Curve("Sideset 1", 1) = {1};
//+
Physical Curve("Sideset 2", 2) = {2};
//+
Physical Curve("Sideset 3", 3) = {3};
//+
Physical Curve("Sideset 4", 4) = {4};
