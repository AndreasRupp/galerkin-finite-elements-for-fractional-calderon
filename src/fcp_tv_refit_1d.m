function [q,info] = fcp_tv_refit_1d(model,sol,alpha_q,q_tv,threshold,bounds)
if nargin<5, threshold=0.5; end
if nargin<6, bounds=[0,1]; end
validateattributes(threshold,{'numeric'},{'scalar','real','finite'});
validateattributes(bounds,{'numeric'},{'vector','numel',2,'real','finite'});
if bounds(1)>bounds(2), error('fcp:TVBounds','The coefficient bounds are not ordered.'); end
validateattributes(alpha_q,{'numeric'},{'scalar','positive','finite'});
n=numel(model.qarea); q_tv=q_tv(:);
if model.mesh.dim~=1 || model.qmesh.dim~=1 || numel(q_tv)~=n || any(~isfinite(q_tv))
    error('fcp:TVRefit','The refit requires one finite TV coefficient per interval.');
end
left=model.qmesh.nodes(model.qmesh.elements(:,1),1);
right=model.qmesh.nodes(model.qmesh.elements(:,2),1);
if any(right<=left) || any(diff(left)<=0) || any(abs(right(1:end-1)-left(2:end))>64*eps*max(1,max(abs([left;right]))))
    error('fcp:TVOrder','The coefficient intervals must be ordered and contiguous.');
end
d=sol.coef_u2(:)+alpha_q*model.qarea(:); b=sol.coef_wu(:);
if numel(d)~=n || numel(b)~=n || any(~isfinite([d;b])) || any(d<=0)
    error('fcp:TVRefit','The quadratic coefficient data must be finite and strictly convex.');
end
q=zeros(n,1); selected=find(q_tv>threshold);
info.threshold=threshold; info.bounds=bounds(:).';
info.support_indices=[]; info.support_interval=[];
info.plateau_unbounded=0; info.plateau=0; info.empty_support=isempty(selected);
if isempty(selected), return; end
indices=(selected(1):selected(end)).'; % The detected support is one interval.
level=-sum(b(indices))/sum(d(indices)); % The constant refit retains the alpha_q penalty.
bounded=min(max(level,bounds(1)),bounds(2));
q(indices)=bounded;
info.support_indices=[indices(1),indices(end)];
info.support_interval=[left(indices(1)),right(indices(end))];
info.plateau_unbounded=level; info.plateau=bounded;
info.scalar_gradient=2*sum(d(indices))*bounded+2*sum(b(indices));
end
