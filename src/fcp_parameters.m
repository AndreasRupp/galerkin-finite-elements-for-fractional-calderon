function cfg = fcp_parameters(example,profile)
if nargin<2, profile='standard'; end
cfg.example=example; cfg.profile=profile;
cfg.coefficient_degree=double(example~=2);

cfg.seed=1234; cfg.noise_generator='portable'; % The fixed seed makes the noise reproducible.
cfg.noise_levels=[1e-7,1e-5,1e-3,1e-1];
cfg.noise_scan=logspace(-10,-1,30).'; cfg.repeats=10;
cfg.mesh_beta=0.45; % This exponent is used only to report the mesh condition.

cfg.order=6; cfg.error_order=10; cfg.stiffness_order=6; % Gauss points per coordinate direction.

cfg.alpha_factor=1; cfg.alpha_power=1.5; cfg.alpha_min=1e-12; % alpha = max(delta^(3/2), 1e-12), with a roundoff safeguard.
cfg.alpha_w_factor=1; cfg.alpha_w_power=2/3; cfg.alpha_w_min=100*eps; % The fractional step uses delta^(2/3).
cfg.alpha_q_factor=0.01; cfg.alpha_q_power=1; cfg.alpha_q_min=1e-12; % The coefficient step uses max(0.01*delta, 1e-12).

cfg.protocols={'independent'}; % The data come from a finer forward mesh.
cfg.primary_protocol='independent';
cfg.run_noise_study=example~=2;
cfg.show_log_fit=example~=2;
cfg.check_quadrature=false;
cfg.example2_panels='both';

cfg.tv.exponents=linspace(-2,1,10);
cfg.tv.selection='fixed'; cfg.tv.fixed_index=1; % The TV parameter is 0.01 times the median baseline jump.
cfg.tv.jump_tolerance=0.05; cfg.tv.target_jumps=2; cfg.tv.jump_weight=1000;
cfg.tv.threshold=0.5; cfg.tv.bounds=[0,1]; % The interval refit uses these prior bounds.
cfg.tv.solver=struct('maxit',10000,'abstol',1e-9,'reltol',1e-7,'gaptol',1e-8);

cfg.make_figures=true; cfg.fig_visible='on'; cfg.close_figures=false;
cfg.surface_export='auto';
root=fileparts(fileparts(mfilename('fullpath')));
cfg.output_dir=fullfile(root,'output');
cfg.fig_dir=fullfile(cfg.output_dir,'figures');

if example==1 || example==2
    cfg.dim=1; cfg.s=0.6; cfg.n=800; cfg.nobs=100;
    cfg.prime=sqrt(0.75);
    cfg.datum.centres=[-1.3;1.3]; cfg.datum.radius=0.2; cfg.datum.offset=1;
    cfg.reference.type='graded'; cfg.reference.n=4000;
    cfg.reference.grading=2.5; cfg.reference.order=10; % The reference mesh is graded towards both endpoints.
    if example==1
        cfg.qfun=@(x) 10*max(0.75-x(:,1).^2,0);
        cfg.q_breaks=[-sqrt(0.75);sqrt(0.75)];
    else
        cfg.qfun=@(x) double(abs(x(:,1))<0.5);
        cfg.q_breaks=[-0.5;0.5];
        cfg.prime=0.8; cfg.noise_levels=[1e-8,1e-7,1e-6,1e-5];
    end
elseif example==3
    cfg.dim=2; cfg.s=0.5; cfg.n=80; cfg.nobs=12;
    cfg.prime=0.75;
    cfg.qfun=@(x) 100*max(0.75^2-x(:,1).^2,0).^3.*max(0.75^2-x(:,2).^2,0).^3;
    cfg.q_breaks=[];
    cfg.datum.centres=[-1.6,0;1.6,0;0,-1.6;0,1.6]; cfg.datum.radius=0.45; cfg.datum.offset=0;
    cfg.reference.type='fine'; cfg.reference.n=640; cfg.reference.order=8; % Eight reference cells per reconstruction cell in each direction.
    cfg.noise_levels=1e-8;
    cfg.noise_scan=logspace(-10,-1,20).'; cfg.repeats=15;
    cfg.alpha_factor=0.1;
    cfg.alpha_min=2e-14; cfg.alpha_q_min=1e-13;
else
    error('fcp:example','The example must be 1, 2 or 3.');
end
datum=cfg.datum;
cfg.ffun=@(x) fcp_exterior_datum(x,datum);

if strcmp(profile,'quick')
    cfg.repeats=2; cfg.noise_scan=logspace(-7,-2,6).'; % The quick profile is only a small trial run.
    if cfg.dim==1
        cfg.n=32; cfg.nobs=40; cfg.reference.n=200;
    else
        cfg.n=8; cfg.nobs=6; cfg.reference.n=32;
        cfg.noise_levels=1e-3;
    end
elseif strcmp(profile,'fine')
    cfg.nobs=2*cfg.nobs;
    cfg.order=8; cfg.stiffness_order=8; cfg.repeats=15;
    if cfg.dim==1
        cfg.n=2*cfg.n;
    else
        cfg.n=128; cfg.reference.n=1280;
    end
elseif ~strcmp(profile,'standard')
    error('fcp:profile','The profile must be quick, standard or fine.');
end
end
