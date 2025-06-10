L=14;
s=2;
H=20;
alpha=Pi/6;
C=(H*H-L*L)/2/L;

h=1/4;

Point(1)={-L,0,0,h};
Point(2)={-L+s,0,0,h};
Point(3)={L-s,0,0,h};
Point(4)={L,0,0,h};

Point(5)={C,0,0,h};
Point(6)={-C,0,0,h};
Point(7)={0,H,0,h};
Point(8)={0,Sqrt((C+L-s)*(C+L-s)-C*C),0,h};

Line(1)={1,2};
Circle(2)={2,5,8};
Circle(3)={8,6,3};
Line(4)={3,4};
Circle(5)={4,6,7};
Circle(6)={7,5,1};

Curve Loop(1)={1:6};
Plane Surface(1)={1};

Physical Curve(1)={1,4};
Physical Curve(2)={5,6};
Physical Surface(800)={1};

Delete{Point(5),Point(6)};

Mesh 2;
Save "mesh.m";