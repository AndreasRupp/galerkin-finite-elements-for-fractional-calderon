function value = fcp_triangle_potential(points,vertices)
if size(points,2)~=2 || ~isequal(size(vertices),[3,2])
    error('fcp:triangle','Expected two-dimensional points and three triangle vertices.');
end
if det([vertices(2,:)-vertices(1,:);vertices(3,:)-vertices(1,:)])<=0
    error('fcp:triangle','Triangle vertices must be counterclockwise and nondegenerate.');
end
value=zeros(size(points,1),1);
for edge=1:3
    a=vertices(edge,:); b=vertices(mod(edge,3)+1,:);
    tangent=b-a; length_edge=sqrt(sum(tangent.^2));
    tangent=tangent/length_edge;
    normal=[tangent(2),-tangent(1)];
    displacement=bsxfun(@minus,a,points);
    distance=displacement*normal.';
    r0=sqrt(sum(displacement.^2,2));
    r1=sqrt(sum(bsxfun(@minus,b,points).^2,2));
    active=abs(distance)>64*eps*length_edge;
    denominator=r0(active)+r1(active)-length_edge;
    edge_integral=zeros(size(denominator));
    close=denominator<sqrt(eps)*length_edge;
    edge_integral(~close)=log1p(2*length_edge./denominator(~close));
    if any(close)
        projected=displacement(active,:)*tangent.';
        transverse=abs(distance(active));
        edge_integral(close)=asinh((projected(close)+length_edge)./transverse(close))-asinh(projected(close)./transverse(close));
    end
    value(active)=value(active)+distance(active).*edge_integral;
end
end
