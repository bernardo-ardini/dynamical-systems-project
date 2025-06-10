classdef FunctionSpace < handle
    properties
        geo
        domain
        constraints
        fe

        B
    end

    methods
        function V=FunctionSpace(geometry,fe)
            V.geo=geometry;
            V.constraints=[];
            V.fe=fe;
        end

        function ndof=numberDof(V)
            if V.fe=="P1"
                ndof=V.geo.numvertices;
            elseif V.fe=="P12"
                ndof=V.geo.d*V.geo.numvertices;
            elseif V.fe=="P12b"
                ndof=V.geo.d*(V.geo.numvertices+V.geo.numtriangles);
            end
        end

        function ndof=numberFreeDof(V)
            if V.fe=="P1"
                ndof=V.geo.numvertices-length(V.constrainedVertices());
            elseif V.fe=="P12"
                ndof=V.geo.d*(V.geo.numvertices-length(V.constrainedVertices()));
            elseif V.fe=="P12b"
                ndof=V.geo.d*(V.geo.numvertices+V.geo.numtriangles-length(V.constrainedVertices()));
            end
        end
        
        function vert=constrainedVertices(V)
            vert=V.geo.lines2vertices(V.constraints);
        end

        function edg=nonConstrainedBoundaryEdges(V)
            pcv=sort(V.geo.lines(ismember(V.geo.lines(:,3),V.constraints),1:2),2);
            ce=full(diag(V.geo.vertices2edges(pcv(:,1),pcv(:,2))));
            edg=V.geo.boundaryEdges();
            edg=edg(~ismember(edg,ce));
        end

        function computeBasisFunctions(V)
            if ismember(V.fe,["P1","P12"])
                V.B=cell(V.geo.numtriangles,1);
                for e=1:V.geo.numtriangles
                    V.B{e}=inv([V.geo.vertices(V.geo.triangles(e,:),:)';ones(1,V.geo.d+1)]);
                end
            end
        end
    end
end