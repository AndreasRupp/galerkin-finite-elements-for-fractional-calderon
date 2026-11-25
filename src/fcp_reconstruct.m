function sol = fcp_reconstruct(model,mu,alpha,alpha_w,alpha_q)
if ~isfield(model,'singular_values')
    error('fcp:prepare','The model must be prepared with fcp_inverse_prepare.');
end
validateattributes(alpha,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(alpha_w,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(alpha_q,{'numeric'},{'scalar','real','finite','positive'});
mu=mu(:);
if numel(mu)~=size(model.B,1), error('fcp:data','The data vector has the wrong size.'); end
sig=model.singular_values;
filter=sig./(sig.^2+alpha);
sol.u=model.Rs\(model.V*(filter.*(model.U'*(model.Rw*mu)))); % Tikhonov solution in the weighted singular basis.
rhs=model.A*sol.u+model.bext; % The fractional term includes the known exterior datum.
if isfield(model,'fractional_alpha') && model.fractional_alpha==alpha_w
    Cw=model.fractional_factor;
    sol.w=Cw\(Cw'\rhs);
else
    sol.w=(model.M+alpha_w*model.A)\rhs;
end
uq=model.P*sol.u; wq=model.P*sol.w;
coefficient=fcp_coefficient_solve(model,uq,wq,alpha_q); % Each coefficient cell is solved independently.
sol.q_coefficients=coefficient.coefficients;
sol.q=coefficient.means;
sol.alpha=alpha; sol.alpha_w=alpha_w; sol.alpha_q=alpha_q;
sol.coef_u2=coefficient.u2; sol.coef_wu=coefficient.wu;
lhs=model.B'*(model.W*(model.B*sol.u))+alpha*(model.S*sol.u);
state_rhs=model.B'*(model.W*mu);
sol.residuals.state=norm(lhs-state_rhs)/max(norm(lhs)+norm(state_rhs),realmin);
Kw=model.M*sol.w+alpha_w*(model.A*sol.w);
sol.residuals.fractional=norm(Kw-rhs)/max(norm(Kw)+norm(rhs),realmin);
sol.residuals.coefficient=coefficient.residual;
r=model.B*sol.u-mu;
sol.data_residual=sqrt(max(0,r'*(model.W*r)));
sol.solution_norm=sqrt(max(0,sol.u'*(model.S*sol.u)));
sol.min_solution_prime=min(uq);
end
