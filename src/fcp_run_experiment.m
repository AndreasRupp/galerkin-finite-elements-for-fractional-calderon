function R = fcp_run_experiment(example,profile,panels,options)
if nargin<2 || isempty(profile), profile='standard'; end
if nargin<3 || isempty(panels), panels='both'; end
if nargin<4, options=struct(); end
if ~any(strcmp(panels,{'both','quadratic','tv'})) || (example~=2 && ~strcmp(panels,'both'))
    error('fcp:panels','The quadratic and TV panel choices apply only to Example 4.2.');
end
cfg=apply_options(fcp_parameters(example,profile),options);
cfg.example2_panels=panels;
use_tv=example==2 && ~strcmp(panels,'quadratic');
names=cfg.protocols;
if isempty(names) || ~all(ismember(names,{'matched','independent'})) || ~any(strcmp(names,cfg.primary_protocol))
    error('fcp:protocol','The protocols are matched and independent, and they must contain the primary protocol.');
end
if ~exist(cfg.fig_dir,'dir'), mkdir(cfg.fig_dir); end
fprintf('\nFractional Calderon reconstruction, example %d, profile %s\n',example,profile);
fprintf('Runtime %s\n',version);
fprintf('Dimension %d, s %.6g, interior divisions %d\n',cfg.dim,cfg.s,cfg.n);
if example==2, fprintf('Example 4.2 panels %s\n',panels); end
model=fcp_build_model(cfg,cfg.n);
forward=fcp_forward(model,cfg.qfun);
data=struct();
if any(strcmp(names,'matched'))
    data.matched.mu=forward.mu;
    data.matched.hs=forward.hs;
    data.matched.emodel=error_model(model,fcp_error_prepare(model,forward,cfg));
end
fprintf('DOFs interior %d, observations %d, coefficient cells %d, degree %d\n',size(model.A,1),size(model.B,1),numel(model.qarea),cfg.coefficient_degree);
fprintf('Mesh diameters interior %.6g, observation %.6g, coefficient %.6g\n',model.mesh.h,model.obs.h,model.qmesh.h);
fprintf('Largest singular value of the observation matrix %.6e\n',model.singular_values(1));
fprintf('Forward residual %.3e, SVD reconstruction residual %.3e\n',forward.normal_residual,model.validation.svd_relative_residual);
if cfg.dim==2
    fprintf('Stiffness pair comparison %.3e, unresolved pairs %d\n',model.assembly.maximum_pair_difference,model.assembly.unresolved_pairs);
    fprintf('Exterior load Gauss order %d, consecutive relative difference %.3e\n',model.exterior_load_info.order,model.exterior_load_info.relative_difference);
end
R.mesh=model.mesh; R.qmesh=model.qmesh; R.representative_index=1;
levels=cfg.noise_levels(:);
if cfg.run_noise_study, levels=[levels;cfg.noise_scan(:)]; end
alphas=arrayfun(@(delta) fcp_parameter_rule(cfg,model,delta),levels);
beta=cfg.mesh_beta;
if beta<=0 || beta>=min(cfg.s,0.5)
    error('fcp:meshExponent','The mesh exponent must lie between zero and min(s,1/2).');
end
R.mesh_indicator=struct('beta',beta,'h',model.h,'minimum_alpha',min(alphas), ...
    'maximum_ratio',model.h^beta/sqrt(min(alphas)));
fprintf('Mesh condition indicator: beta %.2f, h %.6e, min alpha %.6e, max h^beta/sqrt(alpha) %.6e\n', ...
    beta,model.h,min(alphas),R.mesh_indicator.maximum_ratio);
fprintf('The mesh inequality in Lemma 3.6 requires a ratio <= 1. The convergence theorem also requires a joint refinement limit.\n');
R.protocols=names; R.primary_protocol=cfg.primary_protocol;
R.mu_norm.matched=ynorm(model,forward.mu);
R.forward.l2_norm=sqrt(max(0,forward.u'*(model.M*forward.u)));
R.forward.hs_norm=hsnorm(model,forward.u);
if any(strcmp(names,'independent'))
    reference=fcp_reference_forward(cfg,model); % The reference is used for the data and error evaluation.
    data.independent.mu=reference.mu;
    data.independent.hs=reference.hs;
    cache=fcp_error_prepare(model,reference,cfg);
    data.independent.emodel=error_model(model,cache);
    R.mu_norm.independent=ynorm(model,reference.mu);
    R.reference=reference.info;
    R.reference.mesh=reference.model.mesh;
    d=forward.u-cache.u_projection;
    R.forward.data_difference=ynorm(model,reference.mu-forward.mu)/R.mu_norm.independent;
    R.forward.l2_error=sqrt(max(0,d'*(cache.M*d)+cache.u_projection_error2));
    R.forward.hs_error=fcp_hs_distance(reference.hs,forward.u);
    fprintf('Reference mesh %s, %d divisions, element sizes %.3e to %.3e, residual %.3e\n',reference.info.type,reference.info.divisions,reference.info.h_min,reference.info.h_max,reference.normal_residual);
    fprintf('Data of the reference solution and of the forward solution on the reconstruction mesh differ by %.6e relative to ||mu_W||_Y = %.6e\n',R.forward.data_difference,R.mu_norm.independent);
    fprintf('Forward solution on the reconstruction mesh versus reference: L2 %.6e (norm %.6e), Hs %.6e (norm %.6e)\n',R.forward.l2_error,reference.info.l2_norm,R.forward.hs_error,reference.info.hs_norm);
    clear reference cache
end
R.true_u=data.(cfg.primary_protocol).emodel.error_cache.u_reference_nodes;

fprintf('\nRepresentative reconstructions (protocol 1 matched data, 2 independent data)\n');
fprintf('protocol relative_noise absolute_noise alpha alpha_w alpha_q L2_u Hs_u L2_q q_max_error state_res frac_res coef_res\n');
for k=1:numel(cfg.noise_levels)
    rel=cfg.noise_levels(k); seed=cfg.seed+k;
    [alpha,aw,aq,rule]=fcp_parameter_rule(cfg,model,rel);
    model=prepare_fractional(model,aw);
    for p=1:numel(names)
        name=names{p}; id=protocol_id(name);
        [mud,absolute]=fcp_add_noise(model,data.(name).mu,rel,seed,cfg.noise_generator);
        sol=fcp_reconstruct(model,mud,alpha,aw,aq);
        sol.delta=rel; sol.absolute_noise=absolute; sol.seed=seed;
        sol.error=fcp_errors(data.(name).emodel,sol,[],cfg);
        err=sol.error;
        if use_tv
            sol=tv_step(sol,model,data.(name).emodel,aq,cfg);
            ti=sol.tv_info; si=ti.solver;
            fprintf('TV raw %d delta %.3e alpha_tv %.3e iterations %d relative_gap %.3e L2_q %.6e\n',id,rel,ti.alpha_tv,si.iterations,si.relative_primal_dual_gap,sol.tv_error.q_l2);
            fprintf('Prior-based refit %d delta %.3e support [%.6g %.6g] level %.6g L2_q %.6e, from quadratic baseline L2_q %.6e\n',id,rel,ti.support_interval(1),ti.support_interval(2),ti.level_final,sol.tv_refit_error.q_l2,sol.baseline_refit_error.q_l2);
            sol.tv_info=rmfield(sol.tv_info,'candidate_solutions');
        end
        runs.(name)(k)=sol;
        row=[id,rel,absolute,alpha,aw,aq,err.u_l2,err.u_hs,err.w_l2,err.q_l2,err.q_linf, ...
            sol.residuals.state,sol.residuals.fractional,sol.residuals.coefficient,sol.data_residual, ...
            double(rule.alpha_lower_bound_active),double(rule.alpha_q_lower_bound_active),seed];
        fprintf('%d %.3e %.3e %.3e %.3e %.3e %.6e %.6e %.6e %.6e %.3e %.3e %.3e\n',row([1:8,10:14]));
    end
end
R.runs=runs.(cfg.primary_protocol);
R.runs_by_protocol=runs;

if cfg.run_noise_study
    nd=numel(cfg.noise_scan); nr=cfg.repeats; np=numel(names);
    for p=1:np
        N.(names{p})=struct('q_l2',zeros(nd,nr),'q_linf',zeros(nd,nr),'u_hs',zeros(nd,nr),'u_l2',zeros(nd,nr), ...
            'total',zeros(nd,1),'noiseless_error',zeros(nd,1),'noise',zeros(nd,1));
    end
    fprintf('\nNoise study, means over %d realizations at each level\n',nr);
    fprintf('protocol relative_noise alpha mean_L2_q std_L2_q mean_q_max_error mean_Hs_u noiseless_Hs_u noise_Hs_u state_bound_ratio fractional_bound_ratio\n');
    for k=1:nd
        rel=cfg.noise_scan(k); [alpha,aw,aq]=fcp_parameter_rule(cfg,model,rel);
        model=prepare_fractional(model,aw);
        for p=1:np
            name=names{p}; id=protocol_id(name);
            clean=fcp_reconstruct(model,data.(name).mu,alpha,aw,aq);
            N.(name).noiseless_error(k)=fcp_hs_distance(data.(name).hs,clean.u);
            state_ratio=0; fractional_ratio=0;
            for j=1:nr
                seed=cfg.seed+1000*k+j; % Each noise level and realization has its own seed.
                [mud,absolute]=fcp_add_noise(model,data.(name).mu,rel,seed,cfg.noise_generator);
                sol=fcp_reconstruct(model,mud,alpha,aw,aq);
                err=fcp_errors(data.(name).emodel,sol,[],cfg);
                du=sol.u-clean.u; dw=sol.w-clean.w; propagated=hsnorm(model,du);
                N.(name).q_l2(k,j)=err.q_l2; N.(name).q_linf(k,j)=err.q_linf;
                N.(name).u_hs(k,j)=err.u_hs; N.(name).u_l2(k,j)=err.u_l2;
                N.(name).total(k)=N.(name).total(k)+err.u_hs/nr;
                N.(name).noise(k)=N.(name).noise(k)+propagated/nr;
                ratio_u=propagated/max(absolute/(2*sqrt(alpha)),realmin);
                ratio_w=sqrt(max(0,dw'*(model.M*dw)))/max(sqrt(max(0,du'*(model.A*du)))/(2*sqrt(aw)),realmin);
                state_ratio=max(state_ratio,ratio_u); fractional_ratio=max(fractional_ratio,ratio_w);
                if ratio_u>1+1e-5 || ratio_w>1+1e-5
                    error('fcp:bound','A regularization perturbation inequality failed.');
                end
            end
            fprintf('%d %.3e %.3e %.6e %.3e %.6e %.6e %.6e %.6e %.3e %.3e\n',id,rel,alpha,mean(N.(name).q_l2(k,:)),std(N.(name).q_l2(k,:)), ...
                mean(N.(name).q_linf(k,:)),N.(name).total(k),N.(name).noiseless_error(k),N.(name).noise(k),state_ratio,fractional_ratio);
        end
    end
    primary=N.(cfg.primary_protocol);
    R.noise.delta=cfg.noise_scan; R.noise.protocol=cfg.primary_protocol;
    R.noise.q_l2=primary.q_l2; R.noise.q_linf=primary.q_linf;
    R.noise.by_protocol=N;
    R.decomposition.delta=cfg.noise_scan; R.decomposition.protocol=cfg.primary_protocol;
    R.decomposition.total=primary.total; R.decomposition.noiseless_error=primary.noiseless_error; R.decomposition.noise=primary.noise;
    if cfg.show_log_fit
        for p=1:np
            name=names{p}; means=mean(N.(name).q_l2,2);
            fit=fcp_log_fit(cfg.noise_scan,means);
            R.noise.fits.(name)=fit;
            fprintf('Empirical fit C|log(delta)|^(-gamma) for %s data on all %d means: C %.6g, gamma %.6g, log-scale R2 %.6g\n',name,numel(means),fit.C,fit.gamma,fit.r2);
        end
        R.noise.log_fit=R.noise.fits.(cfg.primary_protocol);
    end
end
if cfg.check_quadrature
    higher=cfg; higher.order=cfg.order+2; higher.stiffness_order=cfg.stiffness_order+2; higher.prepare=false;
    mh=fcp_build_model(higher,cfg.n);
    R.quadrature.relative_A=norm(model.A-mh.A,'fro')/norm(mh.A,'fro');
    R.quadrature.relative_B=norm(model.B-mh.B,'fro')/norm(mh.B,'fro');
    R.quadrature.relative_load=norm(model.bext-mh.bext)/norm(mh.bext);
    fprintf('\nQuadrature comparison A %.3e B %.3e exterior_load %.3e\n',R.quadrature.relative_A,R.quadrature.relative_B,R.quadrature.relative_load);
end
R.validation=model.validation;
R.figure_files={};
R.plot_config=cfg; % The figures can be regenerated from the results in memory.
R.figure_error='';
if cfg.make_figures
    try
        R.figure_files=fcp_plot_results(R,cfg);
        fprintf('Figures are saved in %s.\n',cfg.fig_dir);
    catch plot_error
        R.figure_error=plot_error.message;
        warning('fcp:figureExport','Example %d was computed, but figure creation failed: %s',example,plot_error.message);
        fprintf('The numerical results and plotting configuration remain available in the returned result.\n');
    end
end
if isempty(R.figure_error)
    fprintf('Completed example %d.\n',example);
else
    fprintf('Completed numerical calculations for example %d. Figure creation remains incomplete.\n',example);
end
end

function sol=tv_step(sol,model,emodel,alpha_q,cfg)
[tv,info]=fcp_tv_reconstruction_1d(model,sol,alpha_q,cfg.tv);
sol.q_tv=tv.raw; sol.q_tv_refit=tv.refit; sol.q_baseline_refit=tv.baseline_refit;
sol.alpha_tv=info.alpha_tv; sol.tv_info=info;
sol.tv_error=coefficient_error(emodel,sol,tv.raw,cfg);
sol.tv_refit_error=coefficient_error(emodel,sol,tv.refit,cfg);
sol.baseline_refit_error=coefficient_error(emodel,sol,tv.baseline_refit,cfg);
count=size(info.candidate_solutions,2); sol.tv_candidate_errors=zeros(count,1);
for j=1:count
    e=coefficient_error(emodel,sol,info.candidate_solutions(:,j),cfg);
    sol.tv_candidate_errors(j)=e.q_l2;
end
end

function err=coefficient_error(emodel,sol,q,cfg)
comparison=sol; comparison.q=q; comparison.q_coefficients=q;
err=fcp_errors(emodel,comparison,[],cfg);
end

function emodel=error_model(model,cache)
emodel.qarea=model.qarea; emodel.error_cache=cache;
end

function id=protocol_id(name)
id=1+double(strcmp(name,'independent'));
end

function value=hsnorm(model,x)
value=sqrt(max(0,x'*(model.S*x)));
end

function value=ynorm(model,x)
value=sqrt(max(0,x'*(model.W*x)));
end

function model=prepare_fractional(model,alpha_w)
model.fractional_alpha=alpha_w;
model.fractional_factor=chol(model.M+alpha_w*model.A);
end

function cfg=apply_options(cfg,options)
fields=fieldnames(options);
for k=1:numel(fields)
    cfg.(fields{k})=options.(fields{k});
end
if isfield(options,'output_dir')
    if ~isfield(options,'fig_dir'), cfg.fig_dir=fullfile(cfg.output_dir,'figures'); end
end
if ~any(strcmp(cfg.protocols,cfg.primary_protocol))
    cfg.primary_protocol=cfg.protocols{1};
end
end
