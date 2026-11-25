function reference = fcp_forward(model,qfun)
order=model.cfg.order;
if model.mesh.dim==2, order=max(order,8); end
cuts=fcp_potential_breaks(model.cfg);
Q=fcp_reaction(model.mesh,qfun,order,cuts);
K=model.A+Q;
[~,flag]=chol(K);
if flag~=0, error('fcp:forward','The forward matrix is not positive definite.'); end
reference.u=-(K\model.bext);
reference.mu=model.B*reference.u;
reference.reaction=Q;
reference.normal_residual=norm(K*reference.u+model.bext)/max(norm(K*reference.u)+norm(model.bext),realmin);
reference.model.mesh=model.mesh;
reference.hs.S=model.S;
reference.hs.Pc=speye(size(model.S,1));
reference.hs.ref=reference.u;
end
