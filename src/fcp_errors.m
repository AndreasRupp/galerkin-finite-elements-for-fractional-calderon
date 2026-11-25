function err = fcp_errors(model,sol,reference,cfg)
if isfield(model,'error_cache')
    cache=model.error_cache;
else
    cache=fcp_error_prepare(model,reference,cfg);
end
ue=sol.u-cache.u_projection; we=sol.w-cache.w_projection;
err.u_l2=sqrt(max(0,ue'*cache.M*ue+cache.u_projection_error2));
err.w_l2=sqrt(max(0,we'*cache.M*we+cache.w_projection_error2));
err.u_hs=NaN;
if isfield(cache,'hs'), err.u_hs=fcp_hs_distance(cache.hs,sol.u); end
if isfield(sol,'q_coefficients')
    coefficients=sol.q_coefficients;
else
    coefficients=sol.q;
end
qe=coefficients-cache.q_projection_coefficients;
err.q_l2=sqrt(max(0,sum(model.qarea.*sum((qe*cache.q_mass_reference).*qe,2))+cache.q_best_error2));
if cache.qdegree==0
    err.q_linf=max(max(abs(sol.q-cache.qlo),abs(sol.q-cache.qhi)));
    err.q_linf_kind='essential maximum';
else
    values=sum(coefficients(cache.qcell,:).*cache.qbasis,2);
    vertex_error=abs(coefficients-cache.q_vertex_true);
    err.q_linf=max([abs(values-cache.qtrue);vertex_error(:)]);
    err.q_linf_kind='sampled maximum';
end
err.q_best_l2=sqrt(max(0,cache.q_best_error2));
err.u_ref_l2=cache.u_ref_l2; err.w_ref_l2=cache.w_ref_l2;
err.q_min=min(coefficients(:)); err.q_max=max(coefficients(:));
err.q_projection=cache.q_projection;
err.u_reference_nodes=cache.u_reference_nodes;
err.reference_forward_residual=cache.reference_forward_residual;
err.reference_mesh_size=cache.reference_mesh_size;
err.error_quadrature_order=cache.error_quadrature_order;
end
