function reference = fcp_reference_fine_2d(cfg,model)
ref=cfg.reference; n=ref.n; nc=model.mesh.n;
if cfg.dim~=2 || abs(cfg.s-0.5)>64*eps
    error('fcp:order','The fine reference is implemented for dimension two and s=1/2.');
end
if mod(n,nc)~=0
    error('fcp:reference','The reference mesh must be a uniform refinement of the reconstruction mesh.');
end
defaults={'tolerance',1e-13;'maxit',4000;'load_points',96;'load_order',6;'order',8};
for k=1:size(defaults,1)
    if ~isfield(ref,defaults{k,1}), ref.(defaults{k,1})=defaults{k,2}; end
end
h=2/n; started=tic;
op=stiffness_operator(n,cfg.stiffness_order); % The fine-grid operator is applied by FFT without a dense matrix.
x=linspace(-1,1,n+1).';
Q=reaction_matrix(cfg.qfun,x,max(ref.order,8));
[b,load_info]=exterior_load(cfg.datum,x,ref.load_points,ref.load_order);
matrix=@(U) op.apply(U)+reshape(Q*U(:),n-1,n-1);
[U,cg]=conjugate_gradients(matrix,-b,op.precondition,ref.tolerance,ref.maxit);
if ~cg.converged
    error('fcp:forward','The conjugate gradient method for the fine reference did not converge.');
end
u=U(:);
AU=op.apply(U); MU=mass_apply(U,h);
reference.u=u;
reference.mu=observation_values(U,x,model.obs.nodes);
reference.normal_residual=cg.relative_residual;
reference.model.mesh=struct('dim',2,'n',n,'h',sqrt(2)*h,'fine',true);
reference.info.type='fine'; reference.info.divisions=n;
reference.info.h_min=sqrt(2)*h; reference.info.h_max=sqrt(2)*h;
reference.info.order=ref.order;
reference.info.assembly='triangle-pair table applied by the fast Fourier transform';
reference.info.iterations=cg.iterations;
reference.info.pair_table=op.info;
reference.info.load=load_info;
reference.info.hs_norm=sqrt(max(0,u'*(AU(:)+MU(:))));
reference.info.l2_norm=sqrt(max(0,u'*MU(:)));

[X,Y]=ndgrid(x(2:n),x(2:n));
P=fcp_evaluate(model.mesh,[X(:),Y(:)]);
Q2=reaction_matrix(@(z) cfg.qfun(z).^2,x,max(ref.order,8));
w_norm2=max(0,u'*(Q2*u)); clear Q2
Mc=model.M; mu_c=P'*MU(:); mw_c=-(P'*(Q*u));
projections.M=Mc;
projections.u_projection=Mc\mu_c;
projections.w_projection=Mc\mw_c;
projections.u_projection_error2=max(0,u'*MU(:)-projections.u_projection'*mu_c);
projections.w_projection_error2=max(0,w_norm2-projections.w_projection'*mw_c);
projections.u_ref_l2=reference.info.l2_norm;
projections.w_ref_l2=sqrt(w_norm2);
step=n/nc; full=zeros(n+1,n+1); full(2:n,2:n)=U;
coarse=full(1:step:end,1:step:end);
projections.u_reference_nodes=coarse(:);
reference.projections=projections;
rhs=P'*(AU(:)+MU(:));
reference.hs.S=model.S;
reference.hs.Pc=speye(size(model.S,1));
if isfield(model,'Rs')
    reference.hs.ref=model.Rs\(model.Rs'\rhs);
else
    reference.hs.ref=model.S\rhs;
end
reference.hs.offset2=max(0,reference.info.hs_norm^2-reference.hs.ref'*rhs);
reference.info.seconds=toc(started);
end

function op=stiffness_operator(n,order)
[T11,T12,op.info]=fcp_pair_table(n,order);
index=mod(-(n-1):(n-1),2*n)+1;
kernel=zeros(2*n,2*n); kernel(index,index)=T11; K11=fft2(kernel);
kernel=zeros(2*n,2*n); kernel(index,index)=T12; K12=fft2(kernel);
clear T11 T12 kernel
h=2/n;
[k1,k2]=ndgrid(1:n-1,1:n-1);
symbol=sqrt(4*sin(pi*k1/(2*n)).^2+4*sin(pi*k2/(2*n)).^2);
op.n=n; op.h=h;
op.apply=@(U) apply_stiffness(U,n,h,K11,K12);
op.precondition=@(R) sine_transform(sine_transform(R)./symbol)*((2/n)^2/h);
end

function Y=apply_stiffness(U,n,h,K11,K12)
F=zeros(n+1,n+1); F(2:n,2:n)=U;
g1=complex(F(2:end,1:end-1)-F(1:end-1,1:end-1),F(2:end,2:end)-F(2:end,1:end-1));
g2=complex(F(2:end,2:end)-F(1:end-1,2:end),F(1:end-1,2:end)-F(1:end-1,1:end-1));
G1=fft2(g1,2*n,2*n); G2=fft2(g2,2*n,2*n);
Z1=ifft2(K11.*G1+K12.*G2); Z2=ifft2(conj(K12).*G1+K11.*G2);
Z1=Z1(1:n,1:n); Z2=Z2(1:n,1:n);
z1x=real(Z1); z1y=imag(Z1); z2x=real(Z2); z2y=imag(Z2);
Y=zeros(n+1,n+1);
Y(2:end,1:end-1)=Y(2:end,1:end-1)+z1x-z1y;
Y(1:end-1,1:end-1)=Y(1:end-1,1:end-1)-z1x-z2y;
Y(2:end,2:end)=Y(2:end,2:end)+z1y+z2x;
Y(1:end-1,2:end)=Y(1:end-1,2:end)-z2x+z2y;
Y=Y(2:n,2:n)*(h/(2*pi));
end

function Y=sine_transform(X)
Y=sine_columns(sine_columns(X).').';
end

function Y=sine_columns(X)
[m,k]=size(X);
E=zeros(2*(m+1),k); E(2:m+1,:)=X; E(m+3:end,:)=-X(m:-1:1,:);
F=fft(E); Y=-imag(F(2:m+1,:))/2;
end

function Y=mass_apply(U,h)
n=size(U,1)+1; F=zeros(n+1,n+1); F(2:n,2:n)=U;
Y=6*F(2:n,2:n)+F(1:n-1,2:n)+F(3:n+1,2:n)+F(2:n,1:n-1)+F(2:n,3:n+1)+F(3:n+1,3:n+1)+F(1:n-1,1:n-1);
Y=Y*(h^2/12);
end

function [lambda,weights]=triangle_rule(order)
[r,wr]=fcp_gauss(order);
[R,T]=ndgrid(r,r); [WR,WT]=ndgrid(wr,wr);
lambda=[1-R(:),R(:).*(1-T(:)),R(:).*T(:)];
weights=2*R(:).*WR(:).*WT(:);
end

function Q=reaction_matrix(qfun,x,order)
n=numel(x)-1; h=x(2)-x(1); area=h^2/2;
[lambda,weights]=triangle_rule(order);
E=zeros(n,n,3,3,2);
for k=1:numel(weights)
    for type=1:2
        if type==1
            ox=lambda(k,2)+lambda(k,3); oy=lambda(k,3);
        else
            ox=lambda(k,2); oy=lambda(k,2)+lambda(k,3);
        end
        [X,Y]=ndgrid(x(1:n)+h*ox,x(1:n)+h*oy);
        q=reshape(qfun([X(:),Y(:)]),n,n)*(area*weights(k));
        for a=1:3
            for b=a:3
                E(:,:,a,b,type)=E(:,:,a,b,type)+q*(lambda(k,a)*lambda(k,b));
            end
        end
    end
end
[I,J]=ndgrid(1:n,1:n);
node=@(di,dj) I(:)+di+(J(:)+dj-1)*(n+1);
vertices={[node(0,0),node(1,0),node(1,1)],[node(0,0),node(1,1),node(0,1)]};
rows=[]; cols=[]; vals=[];
for type=1:2
    for a=1:3
        for b=1:3
            value=E(:,:,min(a,b),max(a,b),type);
            rows=[rows;vertices{type}(:,a)]; cols=[cols;vertices{type}(:,b)]; vals=[vals;value(:)];
        end
    end
end
Q=sparse(rows,cols,vals,(n+1)^2,(n+1)^2);
inner=false(n+1,n+1); inner(2:n,2:n)=true;
Q=Q(inner(:),inner(:)); Q=(Q+Q.')/2;
end

function [b,info]=exterior_load(datum,x,points,order)
n=numel(x)-1; h=x(2)-x(1); area=h^2/2;
[y,wy]=fcp_datum_quadrature(datum,2);
wy=wy.*fcp_exterior_datum(y,datum);
keep=wy~=0; y=y(keep,:); wy=wy(keep);
[t,bw]=chebyshev_points(points);
[TX,TY]=ndgrid(t,t);
C=reshape(kernel_integral([TX(:),TY(:)],y,wy),points,points);
probe=[0.999,0.013;-0.97,0.31;0.2,-0.995;0.41,0.77;-0.6,-0.05;0.985,0.985];
exact=kernel_integral(probe,y,wy); approx=zeros(size(exact));
for k=1:size(probe,1)
    approx(k)=barycentric(probe(k,1),t,bw)*C*barycentric(probe(k,2),t,bw).';
end
info.interpolation_points=points;
info.interpolation_error=max(abs(approx-exact)./abs(exact));
info.quadrature_order=order;
[lambda,weights]=triangle_rule(order);
B=zeros(n+1,n+1);
for k=1:numel(weights)
    for type=1:2
        if type==1
            ox=lambda(k,2)+lambda(k,3); oy=lambda(k,3);
            shifts=[0,0;1,0;1,1];
        else
            ox=lambda(k,2); oy=lambda(k,2)+lambda(k,3);
            shifts=[0,0;1,1;0,1];
        end
        F=(barycentric(x(1:n)+h*ox,t,bw)*C*barycentric(x(1:n)+h*oy,t,bw).')*(area*weights(k));
        for a=1:3
            ix=(1:n)+shifts(a,1); iy=(1:n)+shifts(a,2);
            B(ix,iy)=B(ix,iy)+lambda(k,a)*F;
        end
    end
end
b=-B(2:n,2:n)/(2*pi);
end

function F=kernel_integral(points,y,wy)
F=zeros(size(points,1),1);
for first=1:512:size(points,1)
    ids=first:min(first+511,size(points,1));
    r2=bsxfun(@minus,points(ids,1),y(:,1).').^2+bsxfun(@minus,points(ids,2),y(:,2).').^2;
    if any(r2(:)<=0), error('fcp:separation','The exterior datum must be separated from Omega.'); end
    F(ids)=(r2.^(-3/2))*wy;
end
end

function [t,w]=chebyshev_points(count)
k=(1:count).'; angle=(2*k-1)*pi/(2*count);
t=cos(angle); w=(-1).^(k-1).*sin(angle);
end

function E=barycentric(points,t,w)
points=points(:); difference=bsxfun(@minus,points,t.');
[row,col]=find(difference==0);
difference(difference==0)=1;
E=bsxfun(@rdivide,w.',difference);
E=bsxfun(@rdivide,E,sum(E,2));
for k=1:numel(row)
    E(row(k),:)=0; E(row(k),col(k))=1;
end
end

function [U,info]=conjugate_gradients(matrix,rhs,precondition,tolerance,maxit)
U=zeros(size(rhs)); R=rhs; Z=precondition(R); D=Z;
rz=sum(sum(R.*Z)); norm_rhs=sqrt(sum(sum(rhs.^2)));
info.converged=false; info.iterations=0; info.relative_residual=1;
if norm_rhs==0, info.converged=true; info.relative_residual=0; return; end
for k=1:maxit
    AD=matrix(D);
    step=rz/sum(sum(D.*AD));
    U=U+step*D; R=R-step*AD;
    info.iterations=k;
    if sqrt(sum(sum(R.^2)))<=tolerance*norm_rhs, break; end
    Z=precondition(R); new=sum(sum(R.*Z));
    D=Z+(new/rz)*D; rz=new;
end
R=rhs-matrix(U);
info.relative_residual=sqrt(sum(sum(R.^2)))/norm_rhs;
info.converged=info.relative_residual<=100*tolerance;
end

function mu=observation_values(U,x,points)
n=numel(x)-1; h=x(2)-x(1);
F=zeros(n+1,n+1); F(2:n,2:n)=U;
g1x=(F(2:end,1:end-1)-F(1:end-1,1:end-1))/h; g1y=(F(2:end,2:end)-F(2:end,1:end-1))/h;
g2x=(F(2:end,2:end)-F(1:end-1,2:end))/h; g2y=(F(1:end-1,2:end)-F(1:end-1,1:end-1))/h;
JH=zeros(n,n+1); JH(:,2:end)=g2y; JH(:,1:end-1)=JH(:,1:end-1)-g1y;
JV=zeros(n+1,n); JV(2:end,:)=g1x; JV(1:end-1,:)=JV(1:end-1,:)-g2x;
JD=(-g1x+g1y+g2x-g2y)/sqrt(2);
[X,Y]=ndgrid(x,x);
mu=zeros(size(points,1),1);
for k=1:size(points,1)
    R=hypot(points(k,1)-X,points(k,2)-Y);
    dh=R(1:end-1,:)+R(2:end,:)-h;
    dv=R(:,1:end-1)+R(:,2:end)-h;
    dd=R(1:end-1,1:end-1)+R(2:end,2:end)-sqrt(2)*h;
    if min([min(dh(:)),min(dv(:)),min(dd(:))])<=0
        error('fcp:separation','An observation point lies on a source edge.');
    end
    mu(k)=sum(sum(JH.*log1p(2*h./dh)))+sum(sum(JV.*log1p(2*h./dv)))+sum(sum(JD.*log1p(2*sqrt(2)*h./dd)));
end
mu=mu/(2*pi);
end
