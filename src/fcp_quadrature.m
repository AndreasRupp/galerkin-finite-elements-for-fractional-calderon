function [points,weights,Phi,element,lambda] = fcp_quadrature(mesh,order)
[r,wr]=fcp_gauss(order);
if mesh.dim==1
    lam=[1-r,r]; wref=wr;
else
    [R,T]=ndgrid(r,r); [WR,WT]=ndgrid(wr,wr);
    lam=[1-R(:),R(:).*(1-T(:)),R(:).*T(:)];
    wref=2*R(:).*WR(:).*WT(:);
end
ne=size(mesh.elements,1); nq=size(lam,1); nv=mesh.dim+1;
element=kron((1:ne).',ones(nq,1));
lambda=repmat(lam,ne,1);
weights=kron(mesh.area,wref);
points=zeros(ne*nq,mesh.dim);
rows=[]; cols=[]; vals=[];
for j=1:nv
    node=kron(mesh.elements(:,j),ones(nq,1));
    points=points+bsxfun(@times,lambda(:,j),mesh.nodes(node,:));
    col=mesh.full_to_interior(node); active=col>0;
    rows=[rows;find(active)]; cols=[cols;col(active)]; vals=[vals;lambda(active,j)];
end
Phi=sparse(rows,cols,vals,ne*nq,numel(mesh.interior));
end
