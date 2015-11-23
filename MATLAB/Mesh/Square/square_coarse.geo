// Gmsh project created on Fri Nov 20 11:02:55 2015


lc = 0.5;

// Physical parameters

hw = 1 ;
tw = 2 ; 

// Define points

// Heart
Point(1) = {-hw, -hw, 0, lc} ;
Point(2) = {hw , -hw, 0, lc} ;
Point(3) = {hw ,  hw, 0, lc} ;
Point(4) = {-hw,  hw, 0, lc} ;

// Torso
Point(5) = {-tw, -tw, 0, lc} ;
Point(6) = {tw , -tw, 0, lc} ;
Point(7) = {tw ,  tw, 0, lc} ;
Point(8) = {-tw,  tw, 0, lc} ;


// Define curves

// Heart
Line(1) = {1,2} ;
Line(2) = {2,3} ;
Line(3) = {3,4} ;
Line(4) = {4,1} ;

// Torso
Line(5) = {5,6} ;
Line(6) = {6,7} ;
Line(7) = {7,8} ;
Line(8) = {8,5} ;


// Connect lines

// Heart
Line Loop(9) = {1,2,3,4} ;

// Torso
Line Loop(10) = {5,6,7,8} ;

// Define surfaces

// Heart
Plane Surface(11) = {9} ;

// Torso
Plane Surface(12) = {10,9};

// Define physical entities 
// This is done to flag appropriate regions

// Dirichlet and Neumann boundaries

FLAG_TORSO_BOUNDARY_DIRI = 1 ; 
FLAG_TORSO_BOUNDARY_NEU = 2 ;
FLAG_HEART_BOUNDARY_DIRI = 3 ;

Physical Line(FLAG_HEART_BOUNDARY_DIRI) = {1,2,3,4} ;

Physical Line(FLAG_TORSO_BOUNDARY_DIRI) = {5,6,7} ;

Physical Line(FLAG_TORSO_BOUNDARY_NEU) = {8} ;

// Region surfaces

FLAG_HEART_REGION = 10 ; 
FLAG_TORSO_REGION = 11 ;

Physical Surface(FLAG_HEART_REGION) = {11} ;
Physical Surface(FLAG_TORSO_REGION) = {12} ;
