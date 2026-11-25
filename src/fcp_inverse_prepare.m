function model = fcp_inverse_prepare(model)
model.Rs=chol(model.S);
model.Rw=chol(model.W);
C=(model.Rw*model.B)/model.Rs; % The change of variables incorporates the data and H^s norms.
[U,Ds,V]=svd(C,'econ');
r=min(size(C));
model.U=U(:,1:r); model.V=V(:,1:r); Ds=Ds(1:r,1:r);
model.singular_values=diag(Ds);
model.validation.whitened_condition=model.singular_values(1)/max(model.singular_values(end),realmin);
model.validation.svd_relative_residual=norm(C-model.U*Ds*model.V','fro')/max(norm(C,'fro'),realmin);
end
