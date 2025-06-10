S=0.5;
L=10;
R=7;
r=6;

H=S/5;
h=H;

Point(4)={R,0,0,h};
Point(5)={r,0,0,h};
Point(6)={-r,0,0,h};
Point(7)={-R,0,0,h};

Point(9)={0,0,0,h};
Point(10)={0,R,0,h};
Point(11)={0,r,0,h};

Circle(1)={4,9,10};
Circle(2)={10,9,7};
Line(3)={7,6};
Circle(4)={6,9,11};
Circle(6)={11,9,5};
Line(7)={5,4};

Curve Loop(1)={1:7};
Plane Surface(1)={1};

Physical Curve(1)={3,7};
Physical Curve(2)={1,2};
Physical Surface(800)={1};

Mesh 2;
Save "mesh.m";