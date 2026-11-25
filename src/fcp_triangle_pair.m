function [value,info] = fcp_triangle_pair(K,L,order)
if nargin<3, order=8; end
validateattributes(order,{'numeric'},{'scalar','integer','>=',2});
if ~isequal(size(K),[3,2]) || ~isequal(size(L),[3,2])
    error('fcp:triangle','Each triangle must have three two-dimensional vertices.');
end
if det([K(2,:)-K(1,:);K(3,:)-K(1,:)])<=0 || det([L(2,:)-L(1,:);L(3,:)-L(1,:)])<=0
    error('fcp:triangle','Triangle vertices must be counterclockwise and nondegenerate.');
end
scale=max(triangle_diameter(K),triangle_diameter(L));
near=norm(mean(K,1)-mean(L,1))<=3*scale;
if near
    tolerance=10^(-min(12,max(8,order+3)));
    q=max(8,2*ceil(order/2));
    previous=symmetric_boundary(K,L,q);
    difference=Inf;
    while q<512
        q=min(2*q,512);
        value=symmetric_boundary(K,L,q);
        difference=abs(value-previous)/max(abs(value),realmin);
        if difference<=tolerance, break; end
        previous=value;
    end
    info.order=q;
    info.relative_difference=difference;
    info.converged=difference<=tolerance;
else
    previous=symmetric_area(K,L,order);
    value=symmetric_area(K,L,order+2);
    info.order=order+2;
    info.relative_difference=abs(value-previous)/max(abs(value),realmin);
    info.converged=true;
end
info.near=near;
if ~(isfinite(value) && value>0)
    error('fcp:triangle','The positive triangle-pair integral was not resolved.');
end
end

function value=symmetric_boundary(K,L,q)
value=(boundary_integral(K,L,q)+boundary_integral(L,K,q))/2;
end

function value=boundary_integral(K,L,q)
[t,w]=cached_gauss(q); value=0;
for e=1:3
    a=K(e,:); b=K(mod(e,3)+1,:); edge=b-a;
    ell=norm(edge); tangent=edge/ell; normal=[tangent(2),-tangent(1)];
    points=bsxfun(@plus,a,t*edge);
    for f=1:3
        c=L(f,:); d=L(mod(f,3)+1,:); source=d-c;
        length_source=norm(source); direction=source/length_source;
        normal_source=[direction(2),-direction(1)];
        factor=normal*normal_source.';
        if abs(factor)<8*eps, continue; end
        difference=bsxfun(@minus,c,points);
        z=difference*direction.';
        transverse=abs(difference(:,1)*direction(2)-difference(:,2)*direction(1));
        first=distance_primitive(z,transverse,length_source);
        last=distance_primitive(z+length_source,transverse,length_source);
        value=value-factor*ell*(w.'*(last-first));
    end
end
end

function value=distance_primitive(z,d,scale)
value=0.5*z.*hypot(z,d);
active=d>64*eps*scale;
value(active)=value(active)+0.5*d(active).^2.*asinh(z(active)./d(active));
end

function value=symmetric_area(K,L,q)
value=(area_potential(K,L,q)+area_potential(L,K,q))/2;
end

function value=area_potential(K,L,q)
[t,w]=cached_gauss(q);
[U,V]=ndgrid(t,t); weights=(w*w.').*(1-U);
points=bsxfun(@plus,K(1,:),U(:)*(K(2,:)-K(1,:))+((1-U(:)).*V(:))*(K(3,:)-K(1,:)));
jacobian=abs(det([K(2,:)-K(1,:);K(3,:)-K(1,:)]));
value=jacobian*(weights(:).'*fcp_triangle_potential(points,L));
end

function diameter=triangle_diameter(T)
diameter=max([norm(T(2,:)-T(1,:)),norm(T(3,:)-T(2,:)),norm(T(1,:)-T(3,:))]);
end

function [x,w]=cached_gauss(n)
persistent nodes weights
if isempty(nodes), nodes=cell(512,1); weights=cell(512,1); end
if numel(nodes)<n || isempty(nodes{n})
    [nodes{n},weights{n}]=fcp_gauss(n);
end
x=nodes{n}; w=weights{n};
end
