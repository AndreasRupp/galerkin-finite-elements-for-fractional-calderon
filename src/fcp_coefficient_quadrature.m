function [x,w,element] = fcp_coefficient_quadrature(mesh,qmesh,order)
if mesh.dim==1
    [r,wr]=fcp_gauss(max(2,order));
    x=[]; w=[]; element=[];
    for j=1:size(qmesh.elements,1)
        a=qmesh.nodes(qmesh.elements(j,1)); b=qmesh.nodes(qmesh.elements(j,2));
        inner=mesh.nodes(mesh.nodes>a+64*eps & mesh.nodes<b-64*eps);
        cuts=[a;inner;b];
        for k=1:numel(cuts)-1
            length=cuts(k+1)-cuts(k);
            x=[x;cuts(k)+length*r];
            w=[w;length*wr];
            element=[element;repmat(j,numel(r),1)];
        end
    end
else
    step=2/mesh.n;
    local=(qmesh.nodes+1)/step;
    if max(abs(local(:)-round(local(:))))>1e-9 || abs(qmesh.h-mesh.h)>1e-9
        error('fcp:coefficientMesh','The two-dimensional coefficient mesh must be aligned with the solution triangulation.');
    end
    [x,w,~,element]=fcp_quadrature(qmesh,max(2,order));
end
end
