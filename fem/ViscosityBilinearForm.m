classdef ViscosityBilinearForm < handle
    properties
        functionSpace
        viscosity        
    end

    methods
        function a=ViscosityBilinearForm(V)
            a.functionSpace=V;
        end

        function A=assembleLoc(a,e)
            V=a.functionSpace;
            geo=V.geo;

            assert((V.fe=="P12")||(V.fe=="P12d"));

            if V.fe=="P12"
                X=[ones(3,1),geo.vertices(geo.triangles(e,:),:)];
                M=det(X)*inv(X)';
                Delta=0.5*abs(det(X));
                b=M(:,2);
                c=M(:,3);
                A=1/(4*Delta)*(b*b'+c*c');
                A=A(:);
                A=A*a.viscosity;
            end
        end

        function A=assemble(a)
            V=a.functionSpace;
            geo=V.geo;

            if V.fe=="P12"
                I=zeros(2*3*3,geo.numtriangles);
                J=zeros(2*3*3,geo.numtriangles);
                V=zeros(2*3*3,geo.numtriangles);
    
                for e=1:geo.numtriangles
                    V(:,e)=a.assembleLoc(e);
                    I(:,e)=repmat([geo.triangles(e,1);geo.triangles(e,2);geo.triangles(e,3)],3,1);
                    J(:,e)=[repmat(geo.triangles(e,1),3,1);repmat(geo.triangles(e,2),3,1);repmat(geo.triangles(e,3),3,1)];
                    V(:,geo.numtriangles-1+e)=V(:,e);
                    I(:,geo.numtriangles-1+e)=repmat([geo.numvertices-1+geo.triangles(e,1);geo.numvertices-1+geo.triangles(e,2);geo.numvertices-1+geo.triangles(e,3)],3,1);
                    J(:,geo.numtriangles-1+e)=[repmat(geo.numvertices-1+geo.triangles(e,1),3,1);repmat(geo.numvertices-1+geo.triangles(e,2),3,1);repmat(geo.numvertices-1+geo.triangles(e,3),3,1)];
                end
    
                V=V(:);
                I=I(:);
                J=J(:);
                
                A=sparse(I,J,V);
                A(a.functionSpace.constrainedVertices(),:)=[];
                A(:,a.functionSpace.constrainedVertices())=[];
                A(geo.numvertices-1+a.functionSpace.constrainedVertices(),:)=[];
                A(:,geo.numvertices-1+a.functionSpace.constrainedVertices())=[];
            end
        end
    end
end