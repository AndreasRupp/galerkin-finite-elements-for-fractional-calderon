function reference = fcp_reference_forward(cfg,model)
ref=cfg.reference; s=cfg.s; breaks=fcp_potential_breaks(cfg);
if cfg.dim==2 && strcmp(ref.type,'fine')
    reference=fcp_reference_fine_2d(cfg,model);
    return
end
if cfg.dim==1
    if strcmp(ref.type,'graded')
        mesh=fcp_graded_mesh(ref.n,ref.grading,breaks);
    elseif strcmp(ref.type,'uniform')
        mesh=fcp_mesh(1,ref.n);
    else
        error('fcp:reference','The reference mesh type must be graded or uniform.');
    end
    [A,assembly]=fcp_stiffness_1d(mesh,s);
    reaction_order=ref.order;
else
    if ~strcmp(ref.type,'uniform') || mod(ref.n,model.mesh.n)~=0
        error('fcp:reference','The two-dimensional reference mesh must be a uniform refinement of the reconstruction mesh.');
    end
    if abs(s-0.5)>64*eps
        error('fcp:order','The two-dimensional assembly is implemented for s=1/2.');
    end
    mesh=fcp_mesh(2,ref.n);
    [A,assembly]=fcp_stiffness_2d(mesh,cfg.stiffness_order);
    reaction_order=max(ref.order,8);
end
M=fcp_mass(mesh);
B=fcp_observation(mesh,model.obs,s,ref.order);
b=fcp_exterior_load_exact(mesh,cfg.datum,s,ref.order);
Q=fcp_reaction(mesh,cfg.qfun,reaction_order,breaks);
K=A+Q;
[~,flag]=chol(K);
if flag~=0, error('fcp:forward','The reference forward matrix is not positive definite.'); end
u=-(K\b);
reference.u=u;
reference.mu=B*u;
reference.normal_residual=norm(K*u+b)/max(norm(K*u)+norm(b),realmin);
reference.model.mesh=mesh;
reference.info.type=ref.type;
reference.info.divisions=ref.n;
reference.info.h_min=min(element_sizes(mesh));
reference.info.h_max=mesh.h;
reference.info.order=ref.order;
if isfield(assembly,'method')
    reference.info.assembly=assembly.method;
else
    reference.info.assembly='analytic and adaptive triangle-pair integrals';
end
reference.info.hs_norm=sqrt(max(0,u'*(A*u)+u'*(M*u)));
reference.info.l2_norm=sqrt(max(0,u'*(M*u)));
reference.hs=fcp_hs_prepare(model.mesh,mesh,A,M,u,s);
end

function h=element_sizes(mesh)
if mesh.dim==1
    h=mesh.area;
else
    a=mesh.nodes(mesh.elements(:,1),:); b=mesh.nodes(mesh.elements(:,2),:); c=mesh.nodes(mesh.elements(:,3),:);
    h=max([sqrt(sum((a-b).^2,2)),sqrt(sum((b-c).^2,2)),sqrt(sum((c-a).^2,2))],[],2);
end
end
