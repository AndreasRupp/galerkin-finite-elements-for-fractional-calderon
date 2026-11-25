function mesh = fcp_mesh(dim,n)
validateattributes(n,{'numeric'},{'scalar','integer','>=',2});
if dim == 1
    nodes = linspace(-1,1,n+1).';
    elements = [(1:n).',(2:n+1).'];
elseif dim == 2
    [X,Y] = ndgrid(linspace(-1,1,n+1));
    nodes = [X(:),Y(:)];
    [I,J] = ndgrid(1:n,1:n);
    a = I(:)+(J(:)-1)*(n+1);
    b = a+1; c = a+(n+1); d = c+1;
    elements = [a,b,d; a,d,c];
else
    error('fcp:dimension','Only dimensions one and two are implemented.');
end
mesh = fcp_mesh_data(nodes,elements);
mesh.interior = find(all(abs(nodes)<1-64*eps,2));
mesh.full_to_interior = zeros(size(nodes,1),1);
mesh.full_to_interior(mesh.interior) = (1:numel(mesh.interior)).';
mesh.n = n;
end
