S=1;
L=4;
HG=30;
l=1;

H=S/6;
h=l/3;

Point(1)={-L,-S,0,H};
Point(2)={L,-S,0,H};
Point(3)={L,0,0,H};
Point(4)={l,0,0,h};
Point(5)={l,HG,0,h};
Point(6)={-l,HG,0,h};
Point(7)={-l,0,0,h};
Point(8)={-L,0,0,H};

For k In {1:7}
   Line(k)={k,k+1};
EndFor

Line(8)={8,1};

Curve Loop(1)={1:8};
Plane Surface(1)={1};

Physical Curve(1)={1};
Physical Curve(2)={5};
Physical Surface(800)={1};

Mesh 2;
Save "mesh.m";