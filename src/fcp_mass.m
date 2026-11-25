function M = fcp_mass(mesh)
nv=mesh.dim+1; ne=size(mesh.elements,1); nd=numel(mesh.interior);
rows=zeros(ne*nv*nv,1); cols=rows; vals=rows; at=0;
for i=1:nv
    for j=1:nv
        a=mesh.full_to_interior(mesh.elements(:,i));
        b=mesh.full_to_interior(mesh.elements(:,j));
        take=a>0 & b>0; nt=sum(take); ids=at+(1:nt);
        rows(ids)=a(take); cols(ids)=b(take);
        vals(ids)=mesh.area(take)*(1+(i==j))/(nv*(nv+1)); at=at+nt; % Exact local P1 mass entries.
    end
end
M=sparse(rows(1:at),cols(1:at),vals(1:at),nd,nd);
end
