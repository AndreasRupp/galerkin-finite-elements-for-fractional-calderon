function cache = fcp_error_prepare(model,reference,cfg)
rm=reference.model;
order=max(cfg.error_order,14*(cfg.dim==2));
if isfield(reference,'projections')
    names=fieldnames(reference.projections);
    for k=1:numel(names), cache.(names{k})=reference.projections.(names{k}); end
else
    if cfg.dim==1
        cuts=unique([model.mesh.nodes;rm.mesh.nodes;target_breaks(cfg)]);
        cuts=cuts(cuts>=-1 & cuts<=1);
        geom=fcp_mesh_data(cuts,[(1:numel(cuts)-1)',(2:numel(cuts))']);
        [x,w]=fcp_quadrature(geom,order);
    else
        [x,w]=fcp_quadrature(rm.mesh,order);
    end
    Pc=fcp_evaluate(model.mesh,x); Pr=fcp_evaluate(rm.mesh,x);
    ur=Pr*reference.u; qr=cfg.qfun(x); wr=-qr(:).*ur;
    cache.M=Pc'*spdiags(w,0,numel(w),numel(w))*Pc;
    cache.M=(cache.M+cache.M')/2;
    cache.u_projection=cache.M\(Pc'*(w.*ur));
    cache.w_projection=cache.M\(Pc'*(w.*wr));
    cache.u_projection_error2=sum(w.*(Pc*cache.u_projection-ur).^2);
    cache.w_projection_error2=sum(w.*(Pc*cache.w_projection-wr).^2);
    cache.u_ref_l2=sqrt(max(0,sum(w.*ur.^2)));
    cache.w_ref_l2=sqrt(max(0,sum(w.*wr.^2)));
    cache.u_reference_nodes=fcp_evaluate(rm.mesh,model.mesh.nodes)*reference.u;
end
[xq,wq,cell]=potential_quadrature(model.qmesh,cfg,order);
qt=cfg.qfun(xq); qt=qt(:);
nc=numel(model.qarea);
degree=0;
if isfield(model,'qdegree'), degree=model.qdegree; end
cache.qdegree=degree;
cache.qbasis=fcp_coefficient_basis(model.qmesh,xq,cell,degree);
cache.qcell=cell; cache.qweights=wq; cache.qtrue=qt;
nv=size(cache.qbasis,2);
rhs=zeros(nc,nv);
for j=1:nv
    rhs(:,j)=accumarray(cell,wq.*qt.*cache.qbasis(:,j),[nc,1]);
end
if degree==0
    cache.q_mass_reference=1;
else
    cache.q_mass_reference=(ones(nv)+eye(nv))/(nv*(nv+1));
end
cache.q_projection_coefficients=bsxfun(@rdivide,rhs,model.qarea)/cache.q_mass_reference;
cache.q_projection=mean(cache.q_projection_coefficients,2);
projection=sum(cache.q_projection_coefficients(cell,:).*cache.qbasis,2);
cache.q_best_error2=sum(wq.*(projection-qt).^2);
if degree==0
    [cache.qlo,cache.qhi]=potential_ranges(model.qmesh,cfg);
else
    points=model.qmesh.nodes(model.qmesh.elements(:),:);
    cache.q_vertex_true=reshape(cfg.qfun(points),nc,nv);
end
cache.reference_forward_residual=reference.normal_residual;
cache.reference_mesh_size=rm.mesh.h;
cache.error_quadrature_order=order;
if isfield(reference,'hs'), cache.hs=reference.hs; end
end

function cuts=target_breaks(cfg)
cuts=fcp_potential_breaks(cfg);
end

function [x,w,cell]=potential_quadrature(mesh,cfg,order)
if mesh.dim==2
    [x,w,~,cell]=fcp_quadrature(mesh,order);
    return
end
[r,wr]=fcp_gauss(max(3,order));
x=[]; w=[]; cell=[]; jumps=target_breaks(cfg);
for j=1:size(mesh.elements,1)
    a=mesh.nodes(mesh.elements(j,1)); b=mesh.nodes(mesh.elements(j,2));
    cuts=[a;jumps(jumps>a & jumps<b);b]; cuts=sort(cuts);
    for k=1:numel(cuts)-1
        h=cuts(k+1)-cuts(k);
        x=[x;cuts(k)+h*r]; w=[w;h*wr];
        cell=[cell;repmat(j,numel(r),1)];
    end
end
end

function [lo,hi]=potential_ranges(mesh,cfg)
ne=size(mesh.elements,1); lo=zeros(ne,1); hi=lo;
if cfg.dim==1
    a=mesh.nodes(mesh.elements(:,1)); b=mesh.nodes(mesh.elements(:,2));
    if cfg.example==1
        far=max(abs(a),abs(b)); near=min(abs(a),abs(b));
        near(a<=0 & b>=0)=0;
        lo=10*max(0.75-far.^2,0); hi=10*max(0.75-near.^2,0);
    elseif cfg.example==2
        inside=a>=-0.5 & b<=0.5;
        outside=b<=-0.5 | a>=0.5;
        lo(inside)=1; hi(~outside)=1;
    else
        error('fcp:maximum','An exact maximum-error rule is not available for this target.');
    end
    return
end
if cfg.example~=3
    error('fcp:maximum','An exact maximum-error rule is not available for this target.');
end
r2=0.75^2;
for j=1:ne
    triangle=mesh.nodes(mesh.elements(j,:),:);
    candidates=triangle;
    for edge=1:3
        a=triangle(edge,:); b=triangle(mod(edge,3)+1,:); v=b-a;
        px=[r2-a(1)^2,-2*a(1)*v(1),-v(1)^2];
        py=[r2-a(2)^2,-2*a(2)*v(2),-v(2)^2];
        p=conv(px,py); dp=(1:4).*p(2:5);
        if any(dp~=0)
            z=roots(fliplr(dp));
            z=real(z(abs(imag(z))<1e-10 & real(z)>-1e-12 & real(z)<1+1e-12));
            z=z(:);
            if ~isempty(z)
                z=min(max(z,0),1);
                candidates=[candidates;bsxfun(@plus,a,z*v)];
            end
        end
    end
    bary=[triangle(2,:)-triangle(1,:);triangle(3,:)-triangle(1,:)]'\(-triangle(1,:)');
    if all(bary>=-1e-12) && sum(bary)<=1+1e-12
        candidates=[candidates;0,0];
    end
    values=100*max(r2-candidates(:,1).^2,0).^3.*max(r2-candidates(:,2).^2,0).^3;
    lo(j)=min(values); hi(j)=max(values);
end
end
