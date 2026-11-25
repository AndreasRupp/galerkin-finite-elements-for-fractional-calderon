function mesh = fcp_mesh_data(nodes,elements)
dim = size(nodes,2); ne = size(elements,1); nv = dim+1;
mesh.nodes = nodes; mesh.elements = elements; mesh.dim = dim;
mesh.area = zeros(ne,1); mesh.grad = zeros(ne,nv,dim);
mesh.h = 0;
if dim == 1
    lengths = nodes(elements(:,2))-nodes(elements(:,1));
    if any(lengths<=0), error('fcp:mesh','Intervals must be positively oriented.'); end
    mesh.area = lengths;
    mesh.grad(:,1,1) = -1./lengths;
    mesh.grad(:,2,1) = 1./lengths;
    mesh.h = max(lengths);
elseif dim == 2
    a=nodes(elements(:,1),:); b=nodes(elements(:,2),:); c=nodes(elements(:,3),:);
    detJ=(b(:,1)-a(:,1)).*(c(:,2)-a(:,2))-(b(:,2)-a(:,2)).*(c(:,1)-a(:,1));
    if any(detJ<=0), error('fcp:mesh','Triangles must be counterclockwise.'); end
    mesh.area=detJ/2;
    mesh.grad(:,1,1)=(b(:,2)-c(:,2))./detJ;
    mesh.grad(:,1,2)=(c(:,1)-b(:,1))./detJ;
    mesh.grad(:,2,1)=(c(:,2)-a(:,2))./detJ;
    mesh.grad(:,2,2)=(a(:,1)-c(:,1))./detJ;
    mesh.grad(:,3,1)=(a(:,2)-b(:,2))./detJ;
    mesh.grad(:,3,2)=(b(:,1)-a(:,1))./detJ;
    mesh.h=max([sqrt(sum((a-b).^2,2));sqrt(sum((b-c).^2,2));sqrt(sum((c-a).^2,2))]);
else
    error('fcp:dimension','Only interval and triangular meshes are supported.');
end
mesh.interior=(1:size(nodes,1)).';
mesh.full_to_interior=mesh.interior;
end
