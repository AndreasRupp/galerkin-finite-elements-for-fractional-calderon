function out = fcp_coefficient_solve(model,uq,wq,alpha_q)
nc=numel(model.qarea);
element=model.qelem; weights=model.qweights;
if isfield(model,'qbasis')
    basis=model.qbasis;
else
    basis=ones(numel(weights),1);
end
nv=size(basis,2);
out.u2=accumarray(element,weights.*uq.^2,[nc,1]);
out.wu=accumarray(element,weights.*uq.*wq,[nc,1]);
if nv==1
    den=out.u2+alpha_q*model.qarea; % Cellwise quadratic minimization for degree zero.
    out.coefficients=-out.wu./den;
    lhs=den.*out.coefficients; rhs=-out.wu;
else
    matrix=zeros(nc,nv,nv); rhs=zeros(nc,nv);
    for i=1:nv
        rhs(:,i)=-accumarray(element,weights.*uq.*wq.*basis(:,i),[nc,1]);
        for j=1:nv
            matrix(:,i,j)=accumarray(element,weights.*uq.^2.*basis(:,i).*basis(:,j),[nc,1])...
                +alpha_q*model.qarea*(1+(i==j))/(nv*(nv+1));
        end
    end
    out.coefficients=zeros(nc,nv); lhs=out.coefficients;
    for k=1:nc
        local=reshape(matrix(k,:,:),nv,nv);
        out.coefficients(k,:)=(local\rhs(k,:)')';
        lhs(k,:)=(local*out.coefficients(k,:)')';
    end
end
out.means=mean(out.coefficients,2);
out.residual=norm(lhs(:)-rhs(:))/max(norm(lhs(:))+norm(rhs(:)),realmin);
end
