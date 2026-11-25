function Q = fcp_reaction(mesh,qfun,order,cuts)
if nargin<4, cuts=[]; end
if mesh.dim==1 && ~isempty(cuts)
    nodes=unique([mesh.nodes;cuts(:)]);
    nodes=nodes(nodes>=mesh.nodes(1) & nodes<=mesh.nodes(end));
    pieces=fcp_mesh_data(nodes,[(1:numel(nodes)-1)',(2:numel(nodes))']);
    [x,w]=fcp_quadrature(pieces,max(3,order));
    Phi=fcp_evaluate(mesh,x);
else
    [x,w,Phi]=fcp_quadrature(mesh,order);
end
if isa(qfun,'function_handle'), q=qfun(x); else, q=qfun; end
if isscalar(q), q=repmat(q,size(w)); end
q=q(:);
if numel(q)~=numel(w), error('fcp:reaction','The coefficient must be evaluated at every quadrature point.'); end
Q=Phi'*spdiags(w.*q,0,numel(w),numel(w))*Phi;
Q=(Q+Q')/2;
end
