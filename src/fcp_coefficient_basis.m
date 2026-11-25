function basis = fcp_coefficient_basis(mesh,x,element,degree)
if degree==0
    basis=ones(size(x,1),1);
elseif degree==1
    vertices=mesh.elements(element,:);
    a=mesh.nodes(vertices(:,1),:);
    b=mesh.nodes(vertices(:,2),:)-a;
    z=x-a;
    if mesh.dim==1
        t=z./b;
        basis=[1-t,t];
    else
        c=mesh.nodes(vertices(:,3),:)-a;
        determinant=b(:,1).*c(:,2)-b(:,2).*c(:,1);
        t=(z(:,1).*c(:,2)-z(:,2).*c(:,1))./determinant;
        v=(b(:,1).*z(:,2)-b(:,2).*z(:,1))./determinant;
        basis=[1-t-v,t,v];
    end
else
    error('fcp:coefficientDegree','The coefficient degree must be zero or one.');
end
end
