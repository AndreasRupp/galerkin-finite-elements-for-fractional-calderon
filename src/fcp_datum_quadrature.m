function [y,w] = fcp_datum_quadrature(datum,dim)
rho=datum.radius; nc=size(datum.centres,1);
if dim==1
    pieces=64; [t,wt]=fcp_gauss(8);
    edges=linspace(-rho,rho,pieces+1).'; h=2*rho/pieces;
    local=reshape(bsxfun(@plus,edges(1:end-1).',h*t),[],1);
    weights=repmat(h*wt,pieces,1);
    y=zeros(numel(local)*nc,1); w=y;
    for j=1:nc
        ids=(j-1)*numel(local)+(1:numel(local));
        y(ids)=datum.centres(j)+local; w(ids)=weights;
    end
elseif dim==2
    nr=40; nt=96; [t,wt]=fcp_gauss(nr);
    r=rho*t; theta=2*pi*(0:nt-1)/nt;
    [R,T]=ndgrid(r,theta); [WR,~]=ndgrid(rho*wt,theta);
    local=[R(:).*cos(T(:)),R(:).*sin(T(:))];
    weights=WR(:).*R(:)*(2*pi/nt);
    y=zeros(size(local,1)*nc,2); w=zeros(size(local,1)*nc,1);
    for j=1:nc
        ids=(j-1)*size(local,1)+(1:size(local,1));
        y(ids,:)=bsxfun(@plus,datum.centres(j,:),local); w(ids)=weights;
    end
else
    error('fcp:dimension','Only dimensions one and two are implemented.');
end
end
