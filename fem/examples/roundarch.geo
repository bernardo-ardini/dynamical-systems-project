S=0.5;
R=7;
r=5;

H=R*1.618;
h=(R-r)/5;

Point(4)={R,0,0,h};
Point(5)={r,0,0,h};
Point(6)={-r,0,0,h};
Point(7)={-R,0,0,h};

Point(9)={0,0,0,h};
Point(10)={0,R,0,h};
Point(11)={0,r,0,h};

Circle(1)={4,9,10};
Circle(2)={10,9,7};
Circle(6)={6,9,11};
Circle(7)={11,9,5};

Point(14)={r,-H,0,h};
Point(15)={R,-H,0,h};

Line(8)={5,14};
Line(9)={14,15};
Line(10)={15,4};

Point(18)={-r,-H,0,h};
Point(19)={-R,-H,0,h};

Line(5)={18,6};
Line(4)={19,18};
Line(3)={7,19};

Curve Loop(1)={1:10};
Plane Surface(1)={1};

Physical Curve(1)={4,9};
Physical Curve(2)={1,2};
Physical Surface(800)={1};

Mesh 2;
Save "mesh.m";