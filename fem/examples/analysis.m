clear;
close all;

addpath("..");

% geometry

system("gmsh beam.geo");

geo=Geometry();
mesh;
geo.vertices=msh.POS(:,1:2);
geo.triangles=msh.TRIANGLES(:,1:3);
geo.lines=msh.LINES;
geo.initialize();

system("rm mesh.m");

% function space

V=FunctionSpace(geo,"P12");
V.constraints=1;
V.computeBasisFunctions();

% koiter analysis

koi=KoiterAnalysis(V);
[uc,lambdac,u,lambda]=koi.performAnalysis(40,"modeNumber",1,"order",2);
fprintf("\nCritical load: %f\n\n",lambdac);

% plot

figure(1);

func=@(S,F) norm(S,"fro");
style='none';

subplot(1,3,1);
patch("Faces",geo.triangles,"Vertices",geo.vertices,'FaceVertexCData',zeros(geo.numtriangles,1),'FaceColor','flat','LineStyle',style);
pbaspect([1,1,1]);
daspect([1,1,1]);
axis off;

subplot(1,3,2);
P=geo.vertices+reshape(uc.dof,[geo.numvertices,2]);
patch("Faces",geo.triangles,"Vertices",P,'FaceVertexCData',koi.computeTrianglesStress(uc,func),'FaceColor','flat','LineStyle',style);
pbaspect([1,1,1]);
daspect([1,1,1]);
axis off;

subplot(1,3,3);
P=geo.vertices+reshape(u.dof,[geo.numvertices,2]);
patch("Faces",geo.triangles,"Vertices",P,'FaceVertexCData',koi.computeTrianglesStress(u,func),'FaceColor','flat','LineStyle',style);
pbaspect([1,1,1]);
daspect([1,1,1]);
axis off;
