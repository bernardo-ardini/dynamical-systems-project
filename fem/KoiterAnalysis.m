classdef KoiterAnalysis < handle
    properties
        functionSpace % fornire un function space P12 con le condizioni di Dirichlet

        W
        S
        D2W
        D3W
        D4W
    end

    methods
        function koi=KoiterAnalysis(V,varargin)
            koi.functionSpace=V;

            p=inputParser;
            addParameter(p,"constitutiveRelation","KirchhoffStVenant");
            addParameter(p,"lambda",120);
            addParameter(p,"mu",80);
            parse(p,varargin{:});

            if p.Results.constitutiveRelation=="KirchhoffStVenant"      
                Q=tensorprod(eye(2),tensorprod(eye(2),tensorprod(eye(2),eye(2))));
                
                lmb=p.Results.lambda;
                mu=p.Results.mu;
                koi.S=@(F) lmb*(trace(F'*F)/2-1)*F+mu*(F*F'*F-F);
                koi.D2W=@(F) lmb*(tensorprod(F,F)+(trace(F'*F)/2-1)*permute(tensorprod(eye(2),eye(2)),[1,3,2,4]))+mu*(permute(tensorprod(eye(2),F'*F),[1,3,2,4])+permute(tensorprod(F,F),[1,4,3,2])+permute(tensorprod(eye(2),F*F'),[3,1,4,2])-permute(tensorprod(eye(2),eye(2)),[1,3,2,4]));
                koi.D3W=@(F) lmb*(permute(tensorprod(eye(2),tensorprod(eye(2),F)),[1,3,5,6,2,4])+permute(tensorprod(eye(2),tensorprod(eye(2),F)),[1,3,2,4,5,6])+permute(tensorprod(eye(2),tensorprod(eye(2),F)),[5,6,1,3,2,4]))+mu*(permute(tensorprod(eye(2),tensorprod(eye(2),F)),[5,4,1,3,2,6])+permute(tensorprod(eye(2),tensorprod(eye(2),F)),[5,4,1,6,2,3])+permute(tensorprod(eye(2),tensorprod(eye(2),F)),[2,4,1,6,5,3])+permute(tensorprod(eye(2),tensorprod(eye(2),F)),[2,6,1,4,5,3])+permute(tensorprod(eye(2),tensorprod(eye(2),F)),[2,6,5,3,1,4])+permute(tensorprod(eye(2),tensorprod(eye(2),F)),[2,4,5,3,1,6]));
                koi.D4W=@(F) lmb*(permute(Q,[5,7,1,3,6,8,2,4])+permute(Q,[6,8,5,7,1,3,2,4])+permute(Q,[6,8,2,4,1,3,5,7]))+...
                    mu*(permute(Q,[5,8,6,4,1,3,2,7])+permute(Q,[5,4,6,8,1,3,2,7])+permute(Q,[3,6,2,8,1,5,4,7])+...
                    +permute(Q,[3,8,2,7,1,5,4,6])+permute(Q,[2,8,3,6,1,5,4,7])+permute(Q,[2,8,3,7,1,5,4,6]));
            end
        end

        function stress=computeTrianglesStress(koi,u,func)
            V=koi.functionSpace;
            geo=V.geo;

            stress=zeros(geo.numtriangles,1);
            
            for e=1:geo.numtriangles
                Du=u.jacobian{e};
                stress(e)=stress(e)+func(koi.S(eye(2)+Du),eye(2)+Du);
            end
        end

        function stress=computeVerticesStress(koi,u,func)
            V=koi.functionSpace;
            geo=V.geo;

            normalization=zeros(geo.numvertices,1);
            stress=zeros(geo.numvertices,1);
            
            for e=1:geo.numtriangles
                Du=u.jacobian{e};
                
                for a=1:3
                    aa=geo.triangles(e,a);

                    normalization(aa)=normalization(aa)+geo.areas(e);
                    stress(aa)=stress(aa)+geo.areas(e)*func(koi.S(eye(2)+Du),eye(2)+Du);
                end
            end

            stress=stress./normalization;
        end

        function Flambda=assembleFlambda(koi,indices,load)
            V=koi.functionSpace;
            geo=V.geo;
            d=geo.d;

            %edg=V.nonConstrainedBoundaryEdges();
            edg=geo.lines2edges(indices);

            if size(edg,1)==0
                Flambda=zeros(V.numberFreeDof(),1);
                return
            end

            I=zeros(d*(d+1),size(edg,1));
            vals=zeros(d*(d+1),size(edg,1));

            cc=1;

            for b=edg'
                e=geo.edges2triangles(b,2);
                mask=abs(geo.triangles2edges(e,:))==b;
                bb=1:3;
                bb=bb(mask);
                
                nu=sign(geo.triangles2edges(e,bb))*geo.normals(b,:)';

                l=geo.lengths(b);
                v=(geo.vertices(geo.edges(b,1),:)+geo.vertices(geo.edges(b,2),:))/2;
                B=V.B{e};

                if isstring(load) && load=="pressure"
                    val=l*nu*[v,1]*B';
                else
                    val=-l*load*[v,1]*B';
                end
                
                val=val';
                vals(:,cc)=val(:);

                c=1;
                for i=1:d
                    for a=1:d+1
                        I(c,cc)=(i-1)*geo.numvertices+geo.triangles(e,a);
                        c=c+1;
                    end
                end

                cc=cc+1;
            end

            vals=vals(:);
            I=I(:);

            vals=[vals;zeros(2*geo.numvertices,1)];
            I=[I;(1:2*geo.numvertices)'];

            Flambda=accumarray(I,vals);
            Flambda([koi.functionSpace.constrainedVertices();geo.numvertices+koi.functionSpace.constrainedVertices()])=[];
        end

        function Fu=assembleFu(koi,u)
            % matrice associata a v|->Fu(u,l\lambda)v

            geo=koi.functionSpace.geo;
            d=geo.d;

            I=zeros(d^2*(d+1)^2,geo.numtriangles);
            J=zeros(d^2*(d+1)^2,geo.numtriangles);
            vals=zeros(d^2*(d+1)^2,geo.numtriangles);

            for e=1:geo.numtriangles
                Du=u.jacobian{e};

                B=koi.functionSpace.B{e};
                B=B(:,1:d);

                val=geo.areas(e)*permute(tensorprod(tensorprod(koi.D2W(eye(d)+Du),B,2,2),B,3,2),[3,1,4,2]);
                vals(:,e)=val(:);

                c=1;
                for j=1:d
                    for b=1:d+1
                        for i=1:d
                            for a=1:d+1
                                I(c,e)=(i-1)*geo.numvertices+geo.triangles(e,a);
                                J(c,e)=(j-1)*geo.numvertices+geo.triangles(e,b);
                                c=c+1;
                            end
                        end
                    end
                end
            end

            vals=vals(:);
            I=I(:);
            J=J(:);

            Fu=sparse(I,J,vals);
            Fu([koi.functionSpace.constrainedVertices();geo.numvertices+koi.functionSpace.constrainedVertices()],:)=[];
            Fu(:,[koi.functionSpace.constrainedVertices();geo.numvertices+koi.functionSpace.constrainedVertices()])=[];
        end

        function Fuu=assembleFuu(koi,u,w)
            % matrice associata a v|->Fuu(u,l\lambda)wv

            geo=koi.functionSpace.geo;
            d=geo.d;

            I=zeros(d^2*(d+1)^2,geo.numtriangles);
            J=zeros(d^2*(d+1)^2,geo.numtriangles);
            vals=zeros(d^2*(d+1)^2,geo.numtriangles);

            for e=1:geo.numtriangles
                Du=u.jacobian{e};
                Dw=w.jacobian{e};

                B=koi.functionSpace.B{e};
                B=B(:,1:d);

                val=geo.areas(e)*permute(tensorprod(tensorprod(tensorprod(koi.D3W(eye(d)+Du),Dw,[1,2],[1,2]),B,2,2),B,3,2),[3,1,4,2]);
                vals(:,e)=val(:);

                c=1;
                for j=1:d
                    for b=1:d+1
                        for i=1:d
                            for a=1:d+1
                                I(c,e)=(i-1)*geo.numvertices+geo.triangles(e,a);
                                J(c,e)=(j-1)*geo.numvertices+geo.triangles(e,b);
                                c=c+1;
                            end
                        end
                    end
                end
            end

            vals=vals(:);
            I=I(:);
            J=J(:);

            Fuu=sparse(I,J,vals);
            Fuu([koi.functionSpace.constrainedVertices();geo.numvertices+koi.functionSpace.constrainedVertices()],:)=[];
            Fuu(:,[koi.functionSpace.constrainedVertices();geo.numvertices+koi.functionSpace.constrainedVertices()])=[];
        end
    
        function Fuuu=assembleFuuu(koi,u,w,v)
            % matrice associata a v|->Fuu(u,l\lambda)wv

            geo=koi.functionSpace.geo;
            d=geo.d;

            I=zeros(d^2*(d+1)^2,geo.numtriangles);
            J=zeros(d^2*(d+1)^2,geo.numtriangles);
            vals=zeros(d^2*(d+1)^2,geo.numtriangles);

            for e=1:geo.numtriangles
                Du=u.jacobian{e};
                Dw=w.jacobian{e};
                Dv=v.jacobian{e};

                B=koi.functionSpace.B{e};
                B=B(:,1:d);

                val=geo.areas(e)*permute(tensorprod(tensorprod(tensorprod(tensorprod(koi.D4W(eye(d)+Du),Dv,[1,2],[1,2]),Dw,[1,2],[1,2]),B,2,2),B,3,2),[3,1,4,2]);
                vals(:,e)=val(:);

                c=1;
                for j=1:d
                    for b=1:d+1
                        for i=1:d
                            for a=1:d+1
                                I(c,e)=(i-1)*geo.numvertices+geo.triangles(e,a);
                                J(c,e)=(j-1)*geo.numvertices+geo.triangles(e,b);
                                c=c+1;
                            end
                        end
                    end
                end
            end

            vals=vals(:);
            I=I(:);
            J=J(:);

            Fuuu=sparse(I,J,vals);
            Fuuu([koi.functionSpace.constrainedVertices();geo.numvertices+koi.functionSpace.constrainedVertices()],:)=[];
            Fuuu(:,[koi.functionSpace.constrainedVertices();geo.numvertices+koi.functionSpace.constrainedVertices()])=[];
        end

        function [uc,lambdac,u,lambda]=performAnalysis(koi,t,varargin)
            % parse inputs

            p=inputParser;
            addParameter(p,"modeNumber",1);
            addParameter(p,"order",2);
            parse(p,varargin{:});

            % function space

            V=koi.functionSpace;

            % reference configuration

            u0=Function(V);
            u0.computeJacobians();

            % fundamental path

            A=koi.assembleFu(u0);
            b=koi.assembleFlambda(2,"pressure");
            
            tol=1e-9;
            maxit=1e9;
            M=ichol(A);
            U0hat=pcg(A,-b,tol,maxit,M,M');
            
            u0hat=Function(V);
            u0hat.fromFreeDof(U0hat);
            u0hat.computeJacobians();
            
            % critical load
            
            B=koi.assembleFuu(u0,u0hat);
            [VV,LL]=eigs(A,-B,p.Results.modeNumber,'smallestabs');
            Vc=VV(:,end);
            lambdac=LL(end,end);
            
            vc=Function(V);
            vc.fromFreeDof(Vc);
            vc.computeJacobians();
            
            uc=Function(V);
            uc.fromFreeDof(lambdac*U0hat);
            uc.computeJacobians();
            
            % first order
            
            C=koi.assembleFuu(uc,vc);
            lambda0dot=-0.5*(Vc'*C*Vc)/(Vc'*C*U0hat);
            
            % second order
            
            if p.Results.order==2
                D=koi.assembleFuuu(uc,vc,vc);
                E=koi.assembleFu(uc);
                W0ddot=-E\(C*(Vc+2*lambda0dot*U0hat));
                lambda0ddot=-(Vc'*C*W0ddot+lambda0dot*U0hat'*C*W0ddot+1/3*Vc'*D*Vc+lambda0dot*Vc'*D*U0hat+lambda0dot^2*U0hat'*D*U0hat)/(Vc'*C*U0hat);
            elseif p.Results.order==1
                W0ddot=0;
                lambda0ddot=0;
            end

            % bifurcated path
            
            lambda=lambdac+lambda0dot*t+0.5*lambda0ddot*t^2;
            U=lambda*U0hat+Vc*t+0.5*W0ddot*t^2;
            
            u=Function(V);
            u.fromFreeDof(U);
            u.computeJacobians();
        end
    end
end