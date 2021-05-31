//+
SetFactory("OpenCASCADE");


// Cylinder
//+
Cylinder(1) = {0, 0, -1.5131, 0, 0, 3.0262, 0.2294, 2*Pi};
//+
Cylinder(2) = {0, 0, -1.5131, 0, 0, 3.0262, 0.2794, 2*Pi};


// Left head
//+
Sphere(3) = {0, 0, -1.5131, 0.2294, -Pi/2, Pi/2, 2*Pi};
//+
Sphere(4) = {0, 0, -1.5131, 0.2794, -Pi/2, Pi/2, 2*Pi};
//+
BooleanDifference{ Volume{4}; Delete; }{ Volume{3}; Delete; }
//+
BooleanDifference{ Volume{4}; Delete; }{ Volume{2}; }


// Right head
//+
Sphere(5) = {0, 0, 1.5131, 0.2294, -Pi/2, Pi/2, 2*Pi};
//+
Sphere(6) = {0, 0, 1.5131, 0.2794, -Pi/2, Pi/2, 2*Pi};
//+
BooleanDifference{ Volume{6}; Delete; }{ Volume{5}; Delete; }
//+
BooleanDifference{ Volume{6}; Delete; }{ Volume{2}; }


// Remove inner cylinder
//+
BooleanDifference{ Volume{2}; Delete; }{ Volume{1}; Delete; }
//+
BooleanFragments{ Volume{4}; Delete; }{ Volume{2}; Delete; }
//+
BooleanFragments{ Volume{6}; Delete; }{ Volume{2}; Delete; }


// Nozzle
//+
Cylinder(7) = {0, 0, 1.7425, 0, 0, 0.1, 0.05, 2*Pi};
//+
BooleanDifference{ Volume{7}; Delete; }{ Volume{6}; }


// Define physical domains
//+
Physical Volume("Cylinder", 1) = {2};
//+
Physical Volume("Head_left", 2) = {4};
//+
Physical Volume("Head_right", 3) = {6};
//+
Physical Volume("Nozzle", 4) = {7};
//+
Physical Surface("Inner", 1) = {20, 15, 18};
//+
Physical Surface("Fixed", 2) = {23};
