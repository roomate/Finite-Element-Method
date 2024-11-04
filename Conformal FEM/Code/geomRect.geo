Mesh.MshFileVersion = 2.2
// definition du pas du maillage
h = 0.04;

// d�finition des points (en 3D, raison pour laquelle il y a un 0 en z)
Point(1) = {0, 0, 0, h};
Point(2) = {4.5, 0, 0, h};
Point(3) = {9, 0, 0, h};
Point(4) = {0, 2, 0, h};
Point(5) = {4.5, 2, 0, h};
Point(6) = {9, 2, 0, h};
// d�finition des segments qui relient les points
Line(1) = {4, 1};
Line(2) = {2, 5};
Line(3) = {6, 3};
Line(4) = {1, 2};
Line(5) = {3, 2};
Line(6) = {5, 4};
Line(7) = {5, 6};
// d�finition des contours ferm�s
Line Loop(1) = {1,4,2,6};
Line Loop(2) = {2,5,3,7};
// d�finition des surfaces � partir contours ferm�s
Plane Surface(1) = {1};
Plane Surface(2) = {2};
// d�finition des �l�ments physiques : pour ces �l�ments, nous pourrons r�cup�rer
//									   les r�f�rences 
Physical Point(1) = {1,2,4,5};
Physical Point(2) = {2,3,5,6};
Physical Line(1) = {4,6,1};
Physical Line(2) = {5,3,7};
Physical Surface(1) = {1};
Physical Surface(2) = {2};
