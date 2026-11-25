function [tv,info] = fcp_tv_reconstruction_1d(model,sol,alpha_q,opts)
if nargin<4, opts=struct(); end
if ~isfield(opts,'exponents'), opts.exponents=linspace(-2,1,10); end
if ~isfield(opts,'selection'), opts.selection='fixed'; end
if ~isfield(opts,'fixed_index'), opts.fixed_index=1; end
if ~isfield(opts,'jump_tolerance'), opts.jump_tolerance=0.05; end
if ~isfield(opts,'target_jumps'), opts.target_jumps=2; end
if ~isfield(opts,'jump_weight'), opts.jump_weight=1000; end
if ~isfield(opts,'threshold'), opts.threshold=0.5; end
if ~isfield(opts,'bounds'), opts.bounds=[0,1]; end
if ~isfield(opts,'solver'), opts.solver=struct(); end
validateattributes(opts.exponents,{'numeric'},{'vector','nonempty','real','finite'});
validateattributes(opts.jump_tolerance,{'numeric'},{'scalar','nonnegative','finite'});
validateattributes(opts.target_jumps,{'numeric'},{'scalar','integer','nonnegative'});
validateattributes(opts.jump_weight,{'numeric'},{'scalar','nonnegative','finite'});
validateattributes(alpha_q,{'numeric'},{'scalar','positive','finite'});
d=sol.coef_u2(:)+alpha_q*model.qarea(:); b=sol.coef_wu(:);
if model.mesh.dim~=1 || isempty(d) || numel(d)~=numel(b) || any(~isfinite([d;b])) || any(d<=0)
    error('fcp:TVData','The TV step requires finite, strictly convex P0 coefficient data.');
end
q_ls=-b./d;
if numel(q_ls)>1, scale=median(abs(diff(q_ls)))+1e-14;
else, scale=1e-14; end
parameters=scale*10.^opts.exponents(:);
if any(~isfinite(parameters)) || any(parameters<0)
    error('fcp:TVParameters','The TV candidates must be finite and nonnegative.');
end
count=numel(parameters);
candidates=zeros(count,9); solutions=zeros(numel(d),count); solvers=cell(count,1);
for k=1:count
    [candidate,solver]=fcp_tv_1d(model,sol,alpha_q,parameters(k),opts.solver);
    if ~solver.converged
        error('fcp:TVConvergence','TV candidate %d did not converge.',k);
    end
    jumps=nnz(abs(diff(candidate))>opts.jump_tolerance);
    distance=norm(candidate-q_ls);
    score=opts.jump_weight*abs(jumps-opts.target_jumps)+distance;
    candidates(k,:)=[parameters(k),jumps,distance,sum(d.*(candidate-q_ls).^2),sum(abs(diff(candidate))), ...
        score,solver.iterations,double(solver.converged),solver.relative_primal_dual_gap];
    solutions(:,k)=candidate; solvers{k}=solver;
end
if strcmp(opts.selection,'fixed')
    validateattributes(opts.fixed_index,{'numeric'},{'scalar','integer','>=',1,'<=',count});
    selected=opts.fixed_index;
elseif strcmp(opts.selection,'two_jump_prior')
    [~,selected]=min(candidates(:,6));
else
    error('fcp:TVSelection','The TV selection rule must be fixed or two_jump_prior.');
end
tv.raw=solutions(:,selected);
[tv.refit,refit]=fcp_tv_refit_1d(model,sol,alpha_q,tv.raw,opts.threshold,opts.bounds);
[tv.baseline_refit,baseline]=fcp_tv_refit_1d(model,sol,alpha_q,q_ls,opts.threshold,opts.bounds);
info.alpha_tv=parameters(selected); info.scale=scale; info.selected_index=selected;
info.selection=opts.selection;
info.candidates=candidates; info.candidate_solutions=solutions;
info.candidate_columns={'alpha_tv','jump_count','baseline_distance_l2','quadratic_functional_increase', ...
    'total_variation','two_jump_score','iterations','converged','relative_primal_dual_gap'};
info.solver=solvers{selected}; info.refit=refit; info.baseline_refit=baseline;
info.options=opts; info.selected_jump_count=candidates(selected,2);
info.support_interval=[NaN,NaN]; info.baseline_support_interval=[NaN,NaN];
if ~refit.empty_support, info.support_interval=refit.support_interval; end
if ~baseline.empty_support, info.baseline_support_interval=baseline.support_interval; end
info.level_unconstrained=refit.plateau_unbounded; info.level_final=refit.plateau;
info.baseline_level_final=baseline.plateau;
end
