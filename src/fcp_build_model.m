function model = fcp_build_model(cfg,n)
model.cfg=cfg;
model.mesh=fcp_mesh(cfg.dim,n);
model.obs=fcp_observation_mesh(cfg.dim,cfg.nobs);
if cfg.dim==1
    [model.A,model.assembly]=fcp_stiffness_1d(model.mesh,cfg.s);
else
    if abs(cfg.s-0.5)>64*eps
        error('fcp:order','The two-dimensional assembly is implemented for s=1/2.');
    end
    [model.A,model.assembly]=fcp_stiffness_2d(model.mesh,cfg.stiffness_order);
end
model.M=fcp_mass(model.mesh);
model.S=model.A+model.M; % The discrete H^s norm includes the L2 term.
[model.B,model.W]=fcp_observation(model.mesh,model.obs,cfg.s,cfg.order); % W is the observation mass matrix.
model.f_nodes=cfg.ffun(model.obs.nodes);
model.f_nodes=model.f_nodes(:);
[model.bext,model.exterior_load_info]=fcp_exterior_load(model.mesh,model.obs,model.f_nodes,cfg.s,cfg.order);
nq=max(2,ceil(cfg.prime*n));
base=fcp_mesh(cfg.dim,nq);
model.qmesh=fcp_mesh_data(cfg.prime*base.nodes,base.elements);
model.qmesh.n=nq;
[model.qx,model.qweights,model.qelem]=fcp_coefficient_quadrature(model.mesh,model.qmesh,cfg.order);
model.P=fcp_evaluate(model.mesh,model.qx);
model.qarea=model.qmesh.area;
model.qdegree=0;
if isfield(cfg,'coefficient_degree'), model.qdegree=cfg.coefficient_degree; end
model.qbasis=fcp_coefficient_basis(model.qmesh,model.qx,model.qelem,model.qdegree);
model.h=max([model.mesh.h,model.obs.h,model.qmesh.h]);
model.validation.symmetry_A=norm(model.A-model.A','fro')/max(norm(model.A,'fro'),realmin);
model.validation.symmetry_M=norm(model.M-model.M','fro')/max(norm(model.M,'fro'),realmin);
model.validation.symmetry_W=norm(model.W-model.W','fro')/max(norm(model.W,'fro'),realmin);
[~,pa]=chol(model.A); [~,ps]=chol(model.S); [~,pw]=chol(model.W);
if pa~=0 || ps~=0 || pw~=0
    error('fcp:positivity','A fractional stiffness or mass matrix failed the positive-definiteness check.');
end
model.validation.cholesky_A=true;
model.validation.cholesky_S=true;
model.validation.cholesky_W=true;
if ~isfield(cfg,'prepare') || cfg.prepare
    model=fcp_inverse_prepare(model);
end
end
