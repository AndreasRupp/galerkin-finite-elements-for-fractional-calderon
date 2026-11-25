function [q,info] = fcp_tv_1d(model,sol,alpha_q,alpha_tv,opts)
if model.mesh.dim~=1, error('fcp:dimension','This TV solver is one-dimensional.'); end
if nargin<5, opts=struct(); end
if ~isfield(opts,'maxit'), opts.maxit=10000; end
if ~isfield(opts,'abstol'), opts.abstol=1e-9; end
if ~isfield(opts,'reltol'), opts.reltol=1e-7; end
if ~isfield(opts,'gaptol'), opts.gaptol=1e-8; end
if ~isfield(opts,'adaptive_rho'), opts.adaptive_rho=true; end
validateattributes(alpha_q,{'numeric'},{'scalar','positive','finite'});
validateattributes(alpha_tv,{'numeric'},{'scalar','nonnegative','finite'});
validateattributes(opts.maxit,{'numeric'},{'scalar','positive','integer','finite'});
validateattributes(opts.abstol,{'numeric'},{'scalar','positive','finite'});
validateattributes(opts.reltol,{'numeric'},{'scalar','positive','finite'});
validateattributes(opts.gaptol,{'numeric'},{'scalar','positive','finite'});
area=model.qarea(:); u2=sol.coef_u2(:); b=sol.coef_wu(:); n=numel(area);
if n==0 || numel(u2)~=n || numel(b)~=n || any(~isfinite([area;u2;b])) ...
        || any(area<=0) || any(u2<0)
    error('fcp:TVData','The TV data must have one finite value per cell and positive cell lengths.');
end
d=u2+alpha_q*area; q=-b./d;
D=sparse([(1:n-1)';(1:n-1)'],[(1:n-1)';(2:n)'],[-ones(n-1,1);ones(n-1,1)],n-1,n);
if n<2 || alpha_tv==0
    info=certificate(d,b,q,D,alpha_tv,zeros(n-1,1));
    info.converged=true; info.iterations=0; info.primal=0; info.dual=0;
    info.rho=0; info.rho_updates=0; info.constant_solution=(n==1);
    return
end
constant=-sum(b)/sum(d);
constant_multiplier=cumsum(2*d*constant+2*b);
constant_multiplier=constant_multiplier(1:end-1);
if alpha_tv>=max(abs(constant_multiplier))
    q=constant*ones(n,1);
    info=certificate(d,b,q,D,alpha_tv,constant_multiplier);
    info.converged=info.relative_primal_dual_gap<=opts.gaptol;
    info.iterations=0; info.primal=0; info.dual=0;
    info.rho=0; info.rho_updates=0; info.constant_solution=true;
    return
end
rho=max(median(2*d),realmin);
if isfield(opts,'rho'), rho=opts.rho; end
validateattributes(rho,{'numeric'},{'scalar','positive','finite'});
DtD=D'*D; H=spdiags(2*d,0,n,n); R=chol(H+rho*DtD);
z=D*q; y=zeros(n-1,1); converged=false; updates=0;
for it=1:opts.maxit
    q=R\(R'\(-2*b+rho*D'*(z-y)));
    Dq=D*q; previous=z; v=Dq+y;
    z=sign(v).*max(abs(v)-alpha_tv/rho,0);
    y=y+Dq-z;
    primal=norm(Dq-z); dual=rho*norm(D'*(z-previous));
    tol_primal=sqrt(n-1)*opts.abstol+opts.reltol*max(norm(Dq),norm(z));
    tol_dual=sqrt(n)*opts.abstol*rho+opts.reltol*norm(rho*D'*y);
    if primal<=tol_primal && dual<=tol_dual
        check=certificate(d,b,q,D,alpha_tv,rho*y);
        if check.relative_primal_dual_gap<=opts.gaptol
            converged=true;
            break
        end
    end
    if opts.adaptive_rho && mod(it,25)==0
        relative_primal=primal/max(tol_primal,realmin);
        relative_dual=dual/max(tol_dual,realmin);
        next_rho=rho;
        if relative_primal>10*relative_dual, next_rho=2*rho;
        elseif relative_dual>10*relative_primal, next_rho=rho/2; end
        if next_rho~=rho
            y=(rho/next_rho)*y; rho=next_rho;
            R=chol(H+rho*DtD); updates=updates+1;
        end
    end
end
info=certificate(d,b,q,D,alpha_tv,rho*y);
info.converged=converged; info.iterations=it;
info.primal=primal; info.dual=dual;
info.rho=rho; info.rho_updates=updates; info.constant_solution=false;
if ~converged
    warning('fcp:TV','The TV solver reached its iteration limit. Its result is marked unconverged.');
end
end

function info=certificate(d,b,q,D,alpha_tv,multiplier)
Dq=D*q;
info.stationarity=norm(2*d.*q+2*b+D'*multiplier)/max(norm(2*d.*q)+norm(2*b)+norm(D'*multiplier),realmin);
info.objective=sum(d.*q.^2+2*b.*q)+alpha_tv*sum(abs(Dq));
info.dual_feasibility=max([0;abs(multiplier)-alpha_tv]);
feasible=min(max(multiplier,-alpha_tv),alpha_tv);
dual_load=2*b+D'*feasible;
info.dual_objective=-sum(dual_load.^2./(4*d));
gap=sum(d.*(q+dual_load./(2*d)).^2)+alpha_tv*sum(abs(Dq))-feasible'*Dq;
scale=sum(d.*q.^2+2*abs(b.*q))+alpha_tv*sum(abs(Dq));
info.primal_dual_gap=max(0,gap);
info.relative_primal_dual_gap=max(0,gap)/max(scale,realmin);
end
