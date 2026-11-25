function files = fcp_plot_original_style(R, cfg)
if ~isfield(cfg,'fig_dir'), cfg.fig_dir=fullfile(pwd,'figures'); end
if ~isfield(cfg,'fig_visible'), cfg.fig_visible='on'; end
if ~isfield(cfg,'close_figures'), cfg.close_figures=false; end
if ~isfield(cfg,'show_log_fit'), cfg.show_log_fit=(cfg.example==1); end
if ~isfield(cfg,'example2_panels'), cfg.example2_panels='both'; end
if ~exist(cfg.fig_dir,'dir'), mkdir(cfg.fig_dir); end
names={'defaultTextInterpreter','defaultAxesTickLabelInterpreter','defaultLegendInterpreter', ...
    'defaultAxesFontSize','defaultLineLineWidth'};
values={'latex','latex','latex',24,1.8};
old=apply_defaults(names,values);
cleanup=onCleanup(@()restore_defaults(names,old));
p=original_parameters(cfg.example);
files={};
if R.mesh.dim==1
    files=plot_1d(R,cfg,p,files);
else
    files=plot_2d(R,cfg,p,files);
end
end

function p=original_parameters(example)
p.fig_width_in=6.9; p.fig_height_in=4.5;
p.font_axes=18; p.font_label=22; p.font_legend=16;
p.axes_line_width=1.0;
p.line_width_true=2.8; p.line_width_rec=2.2; p.line_width_dense=2.7;
p.marker_size=8.8; p.axes_position=[0.12,0.14,0.82,0.78];
p.fig_resolution=400; % Raster figures use 400 dpi.
if example==2
    p.font_label=19; p.font_legend=12;
elseif example==3
    p.fig_width_in=7.2; p.fig_height_in=5.6;
    p.font_axes=23; p.font_label=27; p.font_colorbar=21; p.font_legend=20;
    p.axes_position=[0.16,0.16,0.79,0.77];
end
end

function files=plot_1d(R,cfg,p,files)
x=R.mesh.nodes(:,1); domain=[min(x),max(x)];
runs=R.runs; styles={'-','--',':','-.'};
if cfg.example==1
    [fig,ax]=new_axes(cfg,p,'reconstruction');
    utrue=nodal(R.mesh,R.true_u); allU=utrue;
    h=gobjects(numel(runs)+1,1); labels=cell(numel(runs)+1,1);
    h(1)=plot(ax,x,utrue,'LineWidth',p.line_width_true);
    labels{1}=reference_label(R);
    for k=1:numel(runs)
        u=nodal(R.mesh,runs(k).u); allU=[allU;u];
        h(k+1)=plot(ax,x,u,'LineWidth',p.line_width_rec,'LineStyle',styles{mod(k-1,4)+1});
        labels{k+1}=format_noise_label('u',runs(k));
    end
    xlabel(ax,'$x$','Interpreter','latex','FontUnits','points','FontSize',1.1*p.font_axes);
    ylabel(ax,'$u_{0,h}^{\alpha,\delta}(x)$','Interpreter','latex','FontUnits','points','FontSize',1.1*p.font_axes);
    style_axes(ax,p,domain,nice_ylim(allU),true);
    add_boxed_legend(ax,h,labels,p);
    files=save_plot(fig,cfg,p,'u0_reconstruction.pdf',files);
end

if cfg.example==2
    xt=[domain(1);-0.5;-0.5;0.5;0.5;domain(2)]; qt=[0;0;1;1;0;0];
else
    xt=linspace(domain(1),domain(2),4001).'; qt=true_q(cfg,xt);
end
if cfg.example~=2 || ~strcmp(cfg.example2_panels,'tv')
[fig,ax]=new_axes(cfg,p,'reconstruction');
setappdata(fig,'fcp_coefficient_method','quadratic');
h=gobjects(numel(runs)+1,1); labels=cell(numel(runs)+1,1);
h(1)=plot(ax,xt,qt,'LineWidth',p.line_width_true); labels{1}='$q$';
allQ=qt;
for k=1:numel(runs)
    [xx,qq]=coefficient_curve(R.qmesh,runs(k),domain,'');
    h(k+1)=plot(ax,xx,qq,'LineWidth',p.line_width_rec,'LineStyle',styles{mod(k-1,4)+1});
    labels{k+1}=format_noise_label('q',runs(k)); allQ=[allQ;qq];
end
xlabel(ax,'$x$','Interpreter','latex','FontUnits','points','FontSize',1.1*p.font_axes);
ylabel(ax,'$q_h^{\alpha,\delta}(x)$','Interpreter','latex','FontUnits','points','FontSize',1.1*p.font_axes);
columns=1;
if cfg.example==1, columns=2; end
style_axes(ax,p,domain,nice_ylim(allQ),true);
add_boxed_legend(ax,h,labels,p,columns);
if cfg.example==1, name='q_reconstruction.pdf'; else, name='ex43_q_reconstruction.pdf'; end
files=save_plot(fig,cfg,p,name,files);
end

if cfg.example==2 && ~strcmp(cfg.example2_panels,'quadratic')
    fields={'q_tv','q_tv_refit'}; methods={'tv','tv_refit'};
    outputs={'ex42_tv_q_reconstruction.pdf','ex42_tv_refit_q_reconstruction.pdf'};
    for m=1:2
        [fig,ax]=new_axes(cfg,p,'reconstruction');
        setappdata(fig,'fcp_coefficient_method',methods{m});
        h=gobjects(numel(runs)+1,1); labels=cell(numel(runs)+1,1);
        h(1)=plot(ax,xt,qt,'LineWidth',p.line_width_true); labels{1}='$q$';
        allQ=qt;
        for k=1:numel(runs)
            [xx,qq]=coefficient_curve(R.qmesh,runs(k),domain,fields{m});
            h(k+1)=plot(ax,xx,qq,'LineWidth',p.line_width_rec,'LineStyle',styles{mod(k-1,4)+1});
            labels{k+1}=format_noise_label('q',runs(k)); allQ=[allQ;qq];
        end
        xlabel(ax,'$x$','Interpreter','latex','FontUnits','points','FontSize',1.1*p.font_axes);
        ylabel(ax,'$q_h^{\alpha,\delta}(x)$','Interpreter','latex','FontUnits','points','FontSize',1.1*p.font_axes);
        style_axes(ax,p,domain,nice_ylim(allQ),true);
        add_boxed_legend(ax,h,labels,p);
        files=save_plot(fig,cfg,p,outputs{m},files);
    end
end
if cfg.example==2 || ~isfield(R,'decomposition'), return; end
files=decomposition_plot(R.decomposition,cfg,R.mesh.dim,files);
files=noise_plot(R.noise,cfg,p,'q_stability_loglog.pdf',files);
end

function files=decomposition_plot(D,cfg,dim,files)
p=original_parameters(1);
[fig,ax]=new_axes(cfg,p,'decomposition');
x=D.delta(:); y=[D.total(:),D.noiseless_error(:),D.noise(:)];
h=gobjects(3,1); markers={'o','s','d'}; markerSize=[4.2,4.0,4.0];
markIdx=1:5:numel(x);
for k=1:3
    h(k)=loglog(ax,x,y(:,k),'-','LineWidth',2.4);
end
for k=1:3
    loglog(ax,x(markIdx),y(markIdx,k),markers{k},'MarkerSize',markerSize(k),'LineWidth',0.8, ...
        'Color',get(h(k),'Color'),'MarkerFaceColor',get(h(k),'Color'),'HandleVisibility','off');
end
style_log_axes(ax,p);
set(ax,'TickLength',[0.015 0.015],'Position',[0.13,0.15,0.81,0.78]);
xlim(ax,[min(x),max(x)]); ylim(ax,positive_log_ylim(y,1.45,1e-12));
xlabel(ax,'$\delta$','Interpreter','latex','FontUnits','points','FontSize',1.1*p.font_axes);
if isfield(D,'protocol') && strcmp(D.protocol,'independent')
    target='u_0'; subscript='H^s(\mathrm{I\!R})';
    if dim==2, subscript='H^s(\mathrm{I\!R}^2)'; end
else
    target='u_{0,h}'; subscript='s,R,h';
end
ylabel(ax,['$\|\cdot\|_{',subscript,'}$'],'Interpreter','latex','FontUnits','points','FontSize',1.1*p.font_axes);
labels={ ...
    ['$\|u_{0,h}^{\alpha,\delta}-',target,'\|_{',subscript,'}$'], ...
    ['$\|u_{0,h}^{\alpha,0}-',target,'\|_{',subscript,'}$'], ...
    ['$\|u_{0,h}^{\alpha,\delta}-u_{0,h}^{\alpha,0}\|_{',subscript,'}$']};
add_boxed_legend(ax,h,labels,p,1,18);
if dim==1, name='verify_upgradeA_decomp_ex42.pdf';
else, name='ex43_2d_state_error_decomposition.pdf'; end
files=save_plot(fig,cfg,p,name,files);
end

function files=noise_plot(N,cfg,p,name,files)
[fig,ax]=new_axes(cfg,p,'noise');
x=N.delta(:);
if any(~isfinite(x)) || any(x<=0), error('fcp:plotNoiseLevels','Noise levels must be positive and finite.'); end
normLabel='$\|q_h^{\alpha,\delta}-q\|_{L^2(\Omega^{\prime})}$';
if cfg.example==3, normLabel='$\|q_{h}^{\alpha,\delta}-q\|_{L^2(\Omega^{\prime})}$'; end
series={N.q_l2}; labels={normLabel};
h=gobjects(numel(series),1); yall=[];
styles={'-','-.'}; colors=get(ax,'ColorOrder');
for k=1:numel(series)
    v=series{k};
    if size(v,1)~=numel(x) || any(~isfinite(v(:))) || any(v(:)<0)
        error('fcp:plotNoiseSize','Noise errors must be finite, nonnegative, and have one row per level.');
    end
    y=mean(v,2);
    if any(y<=0), error('fcp:plotZeroError','Logarithmic plots require positive mean errors.'); end
    yall=[yall;y];
    if cfg.example==3, lineWidth=2.2; markerSize=4.5; markerWidth=0.8; stride=2;
    else, lineWidth=p.line_width_dense; markerSize=max(4.0,p.marker_size-3.8); markerWidth=1.0; stride=5; end
    h(k)=loglog(ax,x,y,styles{k},'LineWidth',lineWidth,'Color',colors(k,:));
    idx=1:stride:numel(x);
    loglog(ax,x(idx),y(idx),'o','LineWidth',markerWidth,'MarkerSize',markerSize, ...
        'Color',get(h(k),'Color'),'MarkerFaceColor',get(h(k),'Color'),'HandleVisibility','off');
    if cfg.example==3 && isfield(cfg,'noise_levels')
        marked=false(size(x));
        for j=1:numel(cfg.noise_levels)
            marked=marked | abs(log10(x)-log10(cfg.noise_levels(j)))<1e-12;
        end
        if any(marked)
            loglog(ax,x(marked),y(marked),'o','Color',get(h(k),'Color'), ...
                'MarkerFaceColor','w','MarkerSize',7.0,'LineWidth',1.2,'HandleVisibility','off');
        end
    end
end
if cfg.show_log_fit && isfield(N,'log_fit')
    fit=N.log_fit;
    if isfield(fit,'gamma') && isfinite(fit.gamma) && fit.gamma>0 && isfield(fit,'C') && isfinite(fit.C) && fit.C>0
        Cplot=round(fit.C,1,'significant'); % The curve and legend use the same rounded constant.
        yfit=Cplot*abs(log(x)).^(-fit.gamma);
        fprintf('Dashed logarithmic curve for example %d: rounded C %.6g, fitted gamma %.6g.\n',cfg.example,Cplot,fit.gamma);
        if numel(yfit)~=numel(x) || any(~isfinite(yfit)) || any(yfit<=0)
            error('fcp:plotFitSize','The logarithmic fit must have one positive finite value per noise level.');
        end
        if cfg.example==3, fitWidth=1.6; else, fitWidth=max(1.6,p.line_width_dense-0.8); end
        h(end+1)=loglog(ax,x,yfit,'--','LineWidth',fitWidth,'Color',0.25*[1 1 1]);
        labels{end+1}=sprintf('$C|\\log\\delta|^{-\\gamma}\\quad(\\gamma=%.3g,\\ C=%.3g)$',fit.gamma,Cplot);
        yall=[yall;yfit];
    end
end
if cfg.example==3
    box(ax,'on'); grid(ax,'off');
    set(ax,'XScale','log','YScale','log','XMinorTick','off','YMinorTick','off', ...
        'XMinorGrid','off','YMinorGrid','off','FontSize',p.font_axes,'LineWidth',1.0,'TickLabelInterpreter','latex');
    set(ax,'TickDir','out','Position',p.axes_position);
    ylim(ax,[max(min(yall)/1.25,1e-14),max(yall)*1.15]);
else
    style_log_axes(ax,p); ylim(ax,positive_log_ylim(yall,1.40,1e-12));
end
xlim(ax,[min(x),max(x)]);
ensure_visible_log_ticks(ax);
xlabel(ax,'$\delta$','Interpreter','latex','FontUnits','points','FontSize',p.font_label);
ylabel(ax,normLabel,'Interpreter','latex','FontUnits','points','FontSize',p.font_label);
add_boxed_legend(ax,h,labels,p);
files=save_plot(fig,cfg,p,name,files);
end

function ensure_visible_log_ticks(ax)
drawnow;
limits=get(ax,'YLim'); ticks=get(ax,'YTick');
visible=ticks(ticks>=limits(1) & ticks<=limits(2));
if numel(visible)>=2, return; end
exponents=floor(log10(limits(1))):ceil(log10(limits(2)));
ticks=sort(reshape([1;2;5]*10.^exponents,1,[]));
ticks=ticks(ticks>=limits(1) & ticks<=limits(2));
if numel(ticks)>5, ticks=ticks(round(linspace(1,numel(ticks),5))); end
if numel(ticks)<2
    ticks=linspace(limits(1),limits(2),5); ticks=ticks(2:4);
end
set(ax,'YTick',ticks,'YTickLabel',arrayfun(@(value)sprintf('%.3g',value),ticks,'UniformOutput',false));
end

function files=plot_2d(R,cfg,p,files)
mesh=R.mesh; run=R.runs(1);
mode='auto';
if isfield(cfg,'surface_export'), mode=cfg.surface_export; end
if ~any(strcmp(mode,{'auto','vector','image'}))
    error('fcp:plotExport','The surface export must be auto, vector or image.');
end
cfg.surface_raster=strcmp(mode,'image') || (strcmp(mode,'auto') && isfield(mesh,'n') && mesh.n>40);
if isfield(R,'representative_index'), run=R.runs(R.representative_index); end
qcoef=local_coefficients(R.qmesh,run,'');
[v,f,qrec]=coefficient_surface(R.qmesh,qcoef,4);
qtrue=true_q(cfg,v);
[base,quad,xgrid,ygrid]=outside_grid(mesh,R.qmesh);
qoutside=true_q(cfg,base);
all_values=[0;qtrue;qrec;qoutside];
limits=[floor(2*min(all_values))/2,ceil(2*max(all_values))/2];
if limits(2)==limits(1), limits(2)=limits(1)+0.5; end
for k=1:2
    if k==1
        z=qtrue; base_values=qoutside; qfun=cfg.qfun; label='$q(x,y)$'; name='ex43_2d_a_true_q.pdf';
    else
        z=qrec; base_values=zeros(size(qoutside)); qfun=[]; label='$q_{h}^{\alpha,\delta}(x,y)$'; name='ex43_2d_b_q_rec.pdf';
    end
    [fig,ax]=new_axes(cfg,p,'potential');
    patch('Parent',ax,'Faces',quad,'Vertices',[base,base_values], ...
        'FaceVertexCData',base_values,'FaceColor','interp','EdgeColor',[0.25 0.25 0.25], ...
        'LineWidth',0.25,'Tag','fcp_zero_background');
    if isempty(qfun), coefficient_jump_faces(ax,R.qmesh,qcoef,xgrid,ygrid); end
    patch('Parent',ax,'Faces',f,'Vertices',[v,z],'FaceVertexCData',z, ...
        'FaceColor','interp','EdgeColor','none','Tag','fcp_coefficient_surface');
    coefficient_grid(ax,R.qmesh,qcoef,qfun,xgrid,ygrid);
    cb=surface_style(ax,p,label,limits);
    tickStep=max(1,ceil(diff(limits)/5));
    ticks=tickStep*(ceil(limits(1)/tickStep):floor(limits(2)/tickStep));
    set(ax,'XLim',[xgrid(1),xgrid(end)],'YLim',[ygrid(1),ygrid(end)],'ZLim',limits, ...
        'XTick',[xgrid(1),mean(xgrid([1,end])),xgrid(end)], ...
        'YTick',[ygrid(1),mean(ygrid([1,end])),ygrid(end)],'ZTick',ticks, ...
        'PositionConstraint','innerposition','Units','normalized','SortMethod','depth');
    set(cb,'Ticks',unique([limits,ticks]),'Units','normalized','Position',[0.846,0.111,0.031,0.815]);
    set(ax,'Position',[0.13,0.15,0.66,0.76]);
    files=save_plot(fig,cfg,p,name,files);
end

if isfield(R,'decomposition'), files=decomposition_plot(R.decomposition,cfg,mesh.dim,files); end
ut=nodal(mesh,R.true_u); ur=nodal(mesh,run.u);
umax=max(abs([ut;ur])); if umax==0, umax=1; end
ulimit=nonempty_limits([ut;ur]);
[fig,ax]=new_axes(cfg,p,'surface');
mesh_surface(ax,mesh,ut);
label=reference_label(R);
surface_style(ax,p,[label(1:end-1),'(x,y)$'],[-umax,umax]); zlim(ax,ulimit);
files=save_plot(fig,cfg,p,'ex43_2d_d_u_true.pdf',files);
[fig,ax]=new_axes(cfg,p,'surface');
mesh_surface(ax,mesh,ur);
surface_style(ax,p,'$u_{0,h}^{\alpha,\delta}(x,y)$',[-umax,umax]); zlim(ax,ulimit);
files=save_plot(fig,cfg,p,'ex43_2d_e_u_rec.pdf',files);
if isfield(R,'noise'), files=noise_plot(R.noise,cfg,p,'ex43_2d_g_stability.pdf',files); end
end

function mesh_surface(ax,mesh,values)
nodes=mesh.nodes; n=0;
if isfield(mesh,'n'), n=mesh.n; end
stride=max(1,ceil(n/40));
if stride==1 || size(nodes,1)~=(n+1)^2
    trisurf(mesh.elements,nodes(:,1),nodes(:,2),values,'Parent',ax, ...
        'FaceColor','interp','EdgeColor',[0.25 0.25 0.25],'LineWidth',0.25);
    return
end
trisurf(mesh.elements,nodes(:,1),nodes(:,2),values,'Parent',ax,'FaceColor','interp','EdgeColor','none');
X=reshape(nodes(:,1),n+1,n+1); Y=reshape(nodes(:,2),n+1,n+1); Z=reshape(values,n+1,n+1);
x=[]; y=[]; z=[];
for i=1:stride:n+1
    x=[x;X(i,:).';NaN;X(:,i);NaN]; y=[y;Y(i,:).';NaN;Y(:,i);NaN]; z=[z;Z(i,:).';NaN;Z(:,i);NaN];
end
for k=-(n-stride):stride:(n-stride)
    x=[x;diag(X,k);NaN]; y=[y;diag(Y,k);NaN]; z=[z;diag(Z,k);NaN];
end
plot3(ax,x,y,z,'Color',[0.25 0.25 0.25],'LineWidth',0.25,'HandleVisibility','off');
end

function [fig,ax]=new_axes(cfg,p,kind)
fig=figure('Color','w','Visible',cfg.fig_visible);
apply_fig_size(fig,p);
ax=axes(fig); hold(ax,'on');
set(ax,'FontUnits','points','LabelFontSizeMultiplier',1.1);
set(ax,'ColorOrder',[17 113 190;221 84 0;237 177 32;133 22 209;59 170 50]/255);
setappdata(fig,'fcp_plot_kind',kind);
end

function cb=surface_style(ax,p,zlabelStr,clim)
view(ax,40,28);
axis(ax,'tight'); box(ax,'on'); grid(ax,'on');
colormap(ax,'parula');
set(ax,'FontSize',p.font_axes,'LineWidth',1.0,'TickLabelInterpreter','latex');
xlabel(ax,'$x$','Interpreter','latex','FontUnits','points','FontSize',p.font_label);
ylabel(ax,'$y$','Interpreter','latex','FontUnits','points','FontSize',p.font_label);
zlabel(ax,zlabelStr,'Interpreter','latex','FontUnits','points','FontSize',p.font_label);
set(ax,'TickDir','out');
cb=colorbar(ax); set(cb,'TickLabelInterpreter','latex','FontSize',p.font_colorbar);
caxis(ax,clim);
end

function [x,q]=coefficient_curve(mesh,run,domain,field)
c=local_coefficients(mesh,run,field);
left=mesh.nodes(mesh.elements(:,1),1); right=mesh.nodes(mesh.elements(:,2),1);
[left,order]=sort(left); right=right(order); c=c(order,:);
x=[domain(1);left(1);reshape([left,right].',[],1);right(end);domain(2)]; % Both traces are retained at cell interfaces.
q=[0;0;reshape(c.',[],1);0;0];
end

function c=local_coefficients(mesh,run,field)
ne=size(mesh.elements,1); nv=mesh.dim+1;
if ~isempty(field)
    q=run.(field)(:); c=repmat(q,1,nv);
elseif isfield(run,'q_coefficients') && ~isempty(run.q_coefficients)
    c=run.q_coefficients;
    if isequal(size(c),[ne,1])
        c=repmat(c,1,nv);
    elseif ~isequal(size(c),[ne,nv])
        error('fcp:plotCoefficientSize','Local coefficient values must have one row per cell and either one column or one column per vertex.');
    end
else
    q=run.q(:); c=repmat(q,1,nv);
end
if ~isequal(size(c),[ne,nv]) || any(~isfinite(c(:)))
    error('fcp:plotCoefficientSize','Coefficient values must be finite and match the coefficient mesh.');
end
end

function [v,f,z]=coefficient_surface(mesh,c,n)
ne=size(mesh.elements,1); nr=(n+1)*(n+2)/2;
r=zeros(nr,2); map=zeros(n+1,n+1); k=0;
for i=0:n
    for j=0:n-i
        k=k+1; map(i+1,j+1)=k; r(k,:)=[i,j]/n;
    end
end
local_faces=zeros(n^2,3); k=0;
for i=0:n-1
    for j=0:n-i-1
        k=k+1; local_faces(k,:)=[map(i+1,j+1),map(i+2,j+1),map(i+1,j+2)];
        if i+j<=n-2
            k=k+1; local_faces(k,:)=[map(i+2,j+1),map(i+2,j+2),map(i+1,j+2)];
        end
    end
end
bary=[1-sum(r,2),r];
v=zeros(ne*nr,2); z=zeros(ne*nr,1); f=zeros(ne*n^2,3);
for e=1:ne
    rows=(e-1)*nr+(1:nr);
    v(rows,:)=bary*mesh.nodes(mesh.elements(e,:),:);
    z(rows)=bary*c(e,:).';
    f((e-1)*n^2+(1:n^2),:)=local_faces+(e-1)*nr;
end
end

function coefficient_grid(ax,mesh,c,qfun,xgrid,ygrid)
t=linspace(0,1,5);
edges=[1,2;2,3;3,1];
tol=64*eps(max(1,max(abs(mesh.nodes(:)))));
for k=1:size(edges,1)
    i=edges(k,1); j=edges(k,2);
    a=mesh.nodes(mesh.elements(:,i),:);
    b=mesh.nodes(mesh.elements(:,j),:);
    vertical=abs(a(:,1)-b(:,1))<=tol & ismembertol(a(:,1),xgrid,tol,'DataScale',1);
    horizontal=abs(a(:,2)-b(:,2))<=tol & ismembertol(a(:,2),ygrid,tol,'DataScale',1);
    keep=vertical | horizontal;
    if ~any(keep), continue; end
    x=a(keep,1)*(1-t)+b(keep,1)*t;
    y=a(keep,2)*(1-t)+b(keep,2)*t;
    if isempty(qfun)
        z=c(keep,i)*(1-t)+c(keep,j)*t;
    else
        z=reshape(qfun([x(:),y(:)]),size(x));
    end
    x(:,end+1)=NaN; y(:,end+1)=NaN; z(:,end+1)=NaN;
    plot3(ax,reshape(x.',[],1),reshape(y.',[],1),reshape(z.',[],1), ...
        'Color',[0.25 0.25 0.25],'LineWidth',0.25,'HandleVisibility','off');
end
end

function coefficient_jump_faces(ax,mesh,c,xgrid,ygrid)
[v,f,z]=fcp_surface_jumps(mesh,c,1);
if isempty(f), return; end
patch('Parent',ax,'Faces',f,'Vertices',[v,z],'FaceVertexCData',z, ...
    'FaceColor','interp','EdgeColor','none','Tag','fcp_jump_faces');
pairs=[f(:,[1,2]);f(:,[2,3]);f(:,[3,1])];
a=pairs(:,1); b=pairs(:,2);
tol=64*eps(max(1,max(abs(v(:)))));
vertical=all(abs(v(a,:)-v(b,:))<=tol,2) & abs(z(a)-z(b))>tol;
on_grid=ismembertol(v(a,1),xgrid,tol,'DataScale',1) & ismembertol(v(a,2),ygrid,tol,'DataScale',1);
pairs=unique(sort(pairs(vertical & on_grid,:),2),'rows');
if isempty(pairs), return; end
a=pairs(:,1); b=pairs(:,2); gap=nan(size(a));
x=[v(a,1),v(b,1),gap]; y=[v(a,2),v(b,2),gap]; h=[z(a),z(b),gap];
plot3(ax,reshape(x.',[],1),reshape(y.',[],1),reshape(h.',[],1), ...
    'Color',[0.25 0.25 0.25],'LineWidth',0.25,'HandleVisibility','off');
end

function [vertices,faces,x,y]=outside_grid(mesh,qmesh)
lo=min(qmesh.nodes,[],1); hi=max(qmesh.nodes,[],1);
x=surface_grid(mesh.nodes(:,1),qmesh.nodes(:,1));
y=surface_grid(mesh.nodes(:,2),qmesh.nodes(:,2));
[X,Y]=ndgrid(x,y); [I,J]=ndgrid(1:numel(x)-1,1:numel(y)-1);
a=I(:)+(J(:)-1)*numel(x);
faces=[a,a+1,a+1+numel(x),a+numel(x)];
cx=(x(I(:))+x(I(:)+1))/2; cy=(y(J(:))+y(J(:)+1))/2;
outside=cx<lo(1) | cx>hi(1) | cy<lo(2) | cy>hi(2);
faces=faces(outside,:);
[ids,~,local]=unique(faces(:));
faces=reshape(local,size(faces));
vertices=[X(ids),Y(ids)];
end

function grid=surface_grid(full_nodes,coefficient_nodes)
outer=unique(full_nodes); inner=unique(coefficient_nodes);
stride=max(1,ceil((numel(outer)-1)/40));
outer=unique([outer(1:stride:end);outer(end)]);
inner=unique([inner(1:stride:end);inner(end)]);
tol=64*eps(max(1,max(abs(outer))));
grid=[outer(outer<inner(1)-tol);inner;outer(outer>inner(end)+tol)];
end

function values=true_q(cfg,points)
values=cfg.qfun(points); values=values(:);
end

function label=reference_label(R)
label='$u_{0,h}$';
if isfield(R,'primary_protocol') && strcmp(R.primary_protocol,'independent')
    label='$u_0$';
end
end

function u=nodal(mesh,values)
if numel(values)==size(mesh.nodes,1), u=values(:);
elseif numel(values)==numel(mesh.interior)
    u=zeros(size(mesh.nodes,1),1); u(mesh.interior)=values(:);
else, error('fcp:plotStateSize','State coefficients do not match the mesh.'); end
end

function limits=nonempty_limits(values)
values=values(isfinite(values)); limits=[min(values),max(values)];
if limits(2)<=limits(1), limits(2)=limits(1)+1; end
end

function label=format_noise_label(kind,r)
val=r.delta;
switch lower(kind)
    case 'q'
        label=['$q_h^{\alpha,\delta}\,(\delta=' format_sci_tex(val) ')$'];
    case 'u'
        label=['$u_{0,h}^{\alpha,\delta}\,(\delta=' format_sci_tex(val) ')$'];
    otherwise
        label=['$\delta=' format_sci_tex(val) '$'];
end
end

function s=format_sci_tex(val)
if ~isfinite(val) || val<=0, s='0'; return; end
exponent=round(log10(val));
if abs(log10(val)-exponent)<1e-10
    s=['10^{' num2str(exponent) '}'];
else
    exponent=floor(log10(val));
    s=[sprintf('%.2f',val/10^exponent) '\times 10^{' num2str(exponent) '}'];
end
end

function files=save_plot(fig,cfg,p,name,files)
setappdata(fig,'fcp_output_name',name);
path=fullfile(cfg.fig_dir,name);
content='vector'; resolution=p.fig_resolution;
if isfield(cfg,'surface_raster') && cfg.surface_raster && strcmp(getappdata(fig,'fcp_plot_kind'),'surface')
    content='image'; resolution=400;
end
setappdata(fig,'fcp_export_content',content);
if strcmp(getappdata(fig,'fcp_plot_kind'),'potential')
    set(fig,'PaperUnits','inches','PaperSize',[p.fig_width_in,p.fig_height_in], ...
        'PaperPositionMode','manual','PaperPosition',[0,0,p.fig_width_in,p.fig_height_in]);
    drawnow;
    print(fig,path,'-dpdf','-painters');
else
    save_figure_pdf(fig,path,'ContentType',content,'Resolution',resolution,'Append',false);
end
files{end+1}=path;
if cfg.close_figures, close(fig); end
end

function old=apply_defaults(names,values)
old=cell(size(names));
for k=1:numel(names)
    try
        old{k}=get(groot,names{k});
        set(groot,names{k},values{k});
    catch
        old{k}=[];
    end
end
end

function restore_defaults(names,old)
for k=1:numel(names)
    if ~isempty(old{k})
        try
            set(groot,names{k},old{k});
        catch
        end
    end
end
end

function yl = nice_ylim(data)
data = data(:);
data = data(isfinite(data));
if isempty(data)
    yl = [];
    return;
end
vmin = min(data);
vmax = max(data);
if abs(vmax - vmin) < 1e-14
    delta = 0.05 * max(1, abs(vmax));
    vmin = vmin - delta;
    vmax = vmax + delta;
end
pad = 0.06 * (vmax - vmin);
yl = [vmin - pad, vmax + pad];
if vmin >= -1e-12
    yl(1) = 0;
end
end

function yl = positive_log_ylim(data, highFactor, floorVal)
data = data(:);
data = data(isfinite(data) & data > 0);
if isempty(data)
    yl = [floorVal, 1];
    return;
end
yMin = min(data);
yMax = max(data);
yl = [max(yMin / 1.35, floorVal), yMax * highFactor];
end

function style_axes(ax, p, xLimVec, yLimVec, useStandardXTicks)
grid(ax, 'off');
box(ax, 'on');
set(ax, 'FontSize', p.font_axes, 'LineWidth', p.axes_line_width, 'TickLabelInterpreter', 'latex');
set(ax, 'TickDir', 'out', 'Layer', 'top');
set(ax, 'XMinorTick', 'off', 'YMinorTick', 'off');
set(ax, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
try
    set(ax, 'Position', p.axes_position);
catch
end
if nargin >= 3 && ~isempty(xLimVec)
    xlim(ax, xLimVec);
end
if nargin >= 4 && ~isempty(yLimVec)
    ylim(ax, yLimVec);
end
if nargin >= 5 && useStandardXTicks
    set(ax, 'XTick', [-1, -0.5, 0, 0.5, 1]);
end
end

function style_log_axes(ax, p)
grid(ax, 'off');
box(ax, 'on');
set(ax, 'FontSize', p.font_axes, 'LineWidth', p.axes_line_width, 'TickLabelInterpreter', 'latex');
set(ax, 'TickDir', 'out', 'Layer', 'top');
set(ax, 'XMinorTick', 'off', 'YMinorTick', 'off');
set(ax, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
set(ax, 'XScale', 'log', 'YScale', 'log');
set(ax, 'TickLength', [0.012 0.012]);
try
    set(ax, 'Position', p.axes_position);
catch
end
end

function add_boxed_legend(ax, h, labels, p, columns, font_size)
if nargin<5, columns=1; end
if nargin<6, font_size=p.font_legend; end
valid = isgraphics(h);
h = h(valid);
labels = labels(valid);
if isempty(h)
    return;
end
labels = pad_legend_labels(labels);
lgd = legend(ax, h, labels, 'Interpreter', 'latex', 'Location', 'northeast', 'Orientation', 'vertical', 'Box', 'on', 'FontSize', font_size, 'Color', 'white');
try
    set(lgd, 'EdgeColor', [0.35, 0.35, 0.35], 'LineWidth', 0.75);
catch
end
set(lgd, 'NumColumns', columns);
try
    set(lgd, 'AutoUpdate', 'off');
catch
end
fcp_place_legend(ax, lgd);
end

function labels=pad_legend_labels(labels)
for k=1:numel(labels)
    value=strtrim(labels{k});
    if numel(value)<2 || value(1)~='$' || value(end)~='$'
        error('fcp:legendLabel','Legend labels must be mathematical expressions.');
    end
    labels{k}=['$\quad ',value(2:end-1),' \quad$'];
end
end

function apply_fig_size(fig, p)
set(fig, 'Units', 'inches', 'Position', [1 1 p.fig_width_in p.fig_height_in]);
end

