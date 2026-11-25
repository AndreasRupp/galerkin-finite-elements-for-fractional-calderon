function [B,W] = fcp_observation(mesh,obs,s,order)
if mesh.dim==2
    if abs(s-0.5)>64*eps
        error('fcp:order','The analytic two-dimensional observation is implemented for s=1/2.');
    end
    B=fcp_observation_values_2d(mesh,obs.nodes);
    W=fcp_mass(obs);
    return
end
[y,w,Phi]=fcp_quadrature(mesh,order);
P=spdiags(w,0,numel(w),numel(w))*Phi;
nx=size(obs.nodes,1); B=zeros(nx,numel(mesh.interior)); c=fcp_constant(mesh.dim,s);
for first=1:128:nx
    ids=first:min(first+127,nx); r2=zeros(numel(ids),size(y,1));
    for k=1:mesh.dim
        r2=r2+bsxfun(@minus,obs.nodes(ids,k),y(:,k)').^2;
    end
    if any(r2(:)<=0), error('fcp:separation','Observation and source quadrature points must be separated.'); end
    B(ids,:)=-c*(r2.^(-(mesh.dim+2*s)/2))*P;
end
W=fcp_mass(obs);
end
