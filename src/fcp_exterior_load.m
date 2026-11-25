function [b,info] = fcp_exterior_load(mesh,extmesh,f_nodes,s,order)
if numel(f_nodes)~=numel(extmesh.interior)
    error('fcp:exterior','The exterior nodal vector has the wrong length.');
end
if mesh.dim==2
    if abs(s-0.5)>64*eps
        error('fcp:order','The analytic two-dimensional load is implemented for s=1/2.');
    end
    q=max(6,order);
    if q>=30, error('fcp:quadrature','The initial exterior-load order must be below 30.'); end
    previous=load_2d(mesh,extmesh,f_nodes,q);
    difference=Inf;
    while q<30
        q=q+2;
        b=load_2d(mesh,extmesh,f_nodes,q);
        difference=norm(b-previous)/max(norm(b),realmin);
        if difference<=1e-10, break; end
        previous=b;
    end
    info.order=q; info.relative_difference=difference;
    info.converged=difference<=1e-10;
    info.method='analytic_observation_outer_quadrature';
    if ~info.converged
        error('fcp:quadrature','The separated exterior load did not meet the quadrature comparison by order 30.');
    end
    return
end
info.order=order; info.relative_difference=NaN; info.converged=true;
info.method='fixed_gauss';
[x,wx,Phi]=fcp_quadrature(mesh,order);
[y,wy,Psi]=fcp_quadrature(extmesh,order);
fy=Psi*f_nodes(:); active=abs(fy)>0;
y=y(active,:); wy=wy(active).*fy(active);
load=zeros(size(wx)); c=fcp_constant(mesh.dim,s);
for first=1:128:size(x,1)
    ids=first:min(first+127,size(x,1)); r2=zeros(numel(ids),size(y,1));
    for k=1:mesh.dim
        r2=r2+bsxfun(@minus,x(ids,k),y(:,k)').^2;
    end
    if any(r2(:)<=0), error('fcp:separation','The exterior interpolation must be separated from Omega.'); end
    load(ids)=-c*(r2.^(-(mesh.dim+2*s)/2))*wy;
end
b=Phi'*(wx.*load);
end

function b=load_2d(mesh,extmesh,f_nodes,order)
[y,wy,Psi]=fcp_quadrature(extmesh,order);
fy=Psi*f_nodes(:); active=abs(fy)>0;
y=y(active,:); weights=wy(active).*fy(active);
b=zeros(numel(mesh.interior),1);
for first=1:256:size(y,1)
    ids=first:min(first+255,size(y,1));
    values=fcp_observation_values_2d(mesh,y(ids,:));
    b=b+values'*weights(ids);
end
end
