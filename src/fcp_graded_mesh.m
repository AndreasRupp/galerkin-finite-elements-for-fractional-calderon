function mesh = fcp_graded_mesh(n,grading,breaks)
validateattributes(n,{'numeric'},{'scalar','integer','even','>=',4});
validateattributes(grading,{'numeric'},{'scalar','real','>=',1});
if nargin<3, breaks=[]; end
m=n/2; t=(0:m).'/m; half=1-(1-t).^grading;
nodes=[-flipud(half);half(2:end)];
for k=1:numel(breaks)
    if abs(breaks(k))>=1, continue; end
    [~,id]=min(abs(nodes(2:end-1)-breaks(k)));
    nodes(id+1)=breaks(k);
end
if any(diff(nodes)<=0), error('fcp:mesh','The graded mesh nodes are not increasing.'); end
mesh=fcp_mesh_data(nodes,[(1:n).',(2:n+1).']);
mesh.interior=(2:n).';
mesh.full_to_interior=zeros(n+1,1);
mesh.full_to_interior(mesh.interior)=(1:n-1).';
mesh.n=n;
end
