// Gmsh project created on Fri Nov 20 11:02:55 2015


lc = 0.2;

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

Physical Point(1) = {1,2,3,4} ;
Physical Point(2) = {5,6,7,8} ;


Physical Line("Heart boundary") = {1,2,3,4} ;
Physical Line("Torso external boundary") = {5,6,7,8} ;

Physical Surface("Heart") = {11} ;
Physical Surface("Torso") = {12} ;
