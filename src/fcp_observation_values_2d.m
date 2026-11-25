function B = fcp_observation_values_2d(mesh,points)
if mesh.dim~=2 || size(points,2)~=2
    error('fcp:dimension','A triangular mesh and two-dimensional points are required.');
end
ne=size(mesh.elements,1); nn=size(mesh.nodes,1);
oriented=[mesh.elements(:,[1,2]);mesh.elements(:,[2,3]);mesh.elements(:,[3,1])];
owner=repmat((1:ne).',3,1);
[edges,~,edge_number]=unique(sort(oriented,2),'rows');
boundary_edges=edges(accumarray(edge_number,1)==1,:);
if any(ismember(mesh.interior,boundary_edges(:)))
    error('fcp:boundary','The source basis functions must vanish on the mesh boundary.');
end
vector=mesh.nodes(oriented(:,2),:)-mesh.nodes(oriented(:,1),:);
lengths=sqrt(sum(vector.^2,2)); normal=[vector(:,2),-vector(:,1)]./lengths;
rows=repmat(edge_number,3,1); columns=mesh.elements(owner,:); values=zeros(3*ne,3);
for j=1:3
    values(:,j)=mesh.grad(owner,j,1).*normal(:,1)+mesh.grad(owner,j,2).*normal(:,2);
end
D=sparse(rows,columns(:),values(:),size(edges,1),nn);
D=D(:,mesh.interior);
a=mesh.nodes(edges(:,1),:); b=mesh.nodes(edges(:,2),:); ell=sqrt(sum((b-a).^2,2));
B=zeros(size(points,1),numel(mesh.interior));
block=max(1,min(256,floor(2e6/size(edges,1))));
for first=1:block:size(points,1)
    ids=first:min(first+block-1,size(points,1));
    r0=hypot(bsxfun(@minus,points(ids,1),a(:,1).'),bsxfun(@minus,points(ids,2),a(:,2).'));
    r1=hypot(bsxfun(@minus,points(ids,1),b(:,1).'),bsxfun(@minus,points(ids,2),b(:,2).'));
    denominator=bsxfun(@minus,r0+r1,ell.');
    if any(denominator(:)<=0)
        error('fcp:separation','An observation point lies on a source edge.');
    end
    potential=log1p(bsxfun(@rdivide,2*ell.',denominator));
    B(ids,:)=(potential*D)/(2*pi);
end
end
