function [A,info] = fcp_stiffness_1d(mesh,s)
if mesh.dim~=1, error('fcp:dimension','This assembly requires an interval mesh.'); end
c=fcp_constant(1,s); x=mesh.nodes; e=mesh.elements;
a=x(e(:,1)); b=x(e(:,2)); h=b-a;
if max(abs(h-h(1)))<=128*eps*max(1,max(abs(x))) && size(e,1)==numel(x)-1
    if any(mesh.interior==1 | mesh.interior==numel(x))
        error('fcp:boundary','The fractional stiffness requires zero endpoint values.');
    end
    col=uniform_column(numel(x)-2,h(1),s,c);
    fullA=toeplitz(col);
    ids=mesh.interior-1; A=fullA(ids,ids);
    info.method='exact uniform P1 formula with stable far-field series';
    info.symmetry=0;
    info.exterior_tail='included through the full-space kernel';
    return
end
if abs(s-0.5)<32*eps
    F=@logprimitive;
    K=F(bsxfun(@minus,b,a'))-F(bsxfun(@minus,a,a'))-F(bsxfun(@minus,b,b'))+F(bsxfun(@minus,a,b'));
    scale=-c;
else
    power=3-2*s;
    F=@(z) abs(z).^power/((2-2*s)*(3-2*s));
    K=F(bsxfun(@minus,b,a'))-F(bsxfun(@minus,a,a'))-F(bsxfun(@minus,b,b'))+F(bsxfun(@minus,a,b'));
    scale=c/(2*s*(2*s-1));
end

ne=size(e,1); N=size(x,1);
G=sparse([(1:ne)';(1:ne)'],[e(:,1);e(:,2)],[-1./h;1./h],ne,N);
G=G(:,mesh.interior);
A=scale*(G'*K*G); A=(A+A')/2; % The derivative kernel includes the full exterior tail.
info.method='exact integrated derivative kernel';
info.symmetry=norm(A-A','fro')/max(norm(A,'fro'),realmin);
info.exterior_tail='included through the full-space kernel';
end

function col=uniform_column(n,h,s,c)
p=3-2*s; delta=p-2; col=zeros(n,1);
factor=c*h^(1-2*s)/(2*s*(2-2*s)*(3-2*s));
weights=[1,-4,6,-4,1]; % Fourth differences give the uniform P1 stiffness entries.
for k=0:min(2,n-1)
    z=abs(k+(-2:2)); nz=z>0; term=zeros(size(z));
    if abs(delta)<32*eps
        term(nz)=z(nz).^2.*log(z(nz));
    else
        term(nz)=z(nz).^2.*expm1(delta*log(z(nz)))/delta;
    end
    col(k+1)=factor*sum(weights.*term);
end
for k=3:n-1
    dividedBinomial=p*(p-1)/6; total=0;
    for j=4:160
        dividedBinomial=dividedBinomial*(p-j+1)/j;
        if mod(j,2)==0
            term=2*dividedBinomial*k^p*((2/k)^j-4*(1/k)^j);
            total=total+term;
            if j>=12 && abs(term)<=eps*abs(total), break; end
        end
    end
    col(k+1)=factor*total;
end
end

function y=logprimitive(x)
y=zeros(size(x)); nz=x~=0; z=abs(x(nz));
y(nz)=z.^2.*(0.5*log(z)-0.75);
end
