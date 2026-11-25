function obs = fcp_observation_mesh(dim,n)
validateattributes(n,{'numeric'},{'scalar','integer','positive'});
if dim==1
    xl=linspace(-3,-1.05,n+1).'; xr=linspace(1.05,3,n+1).';
    nodes=[xl;xr]; e=[(1:n).',(2:n+1).'];
    elements=[e;e+n+1];
else
    ngap=max(1,round(2.1*n/1.95));
    a=linspace(-3,-1.05,n+1); b=linspace(-1.05,1.05,ngap+1); c=linspace(1.05,3,n+1);
    axis=[a,b(2:end),c(2:end)].'; N=numel(axis);
    [X,Y]=ndgrid(axis); nodes=[X(:),Y(:)];
    [I,J]=ndgrid(1:N-1,1:N-1);
    centreX=(axis(I(:))+axis(I(:)+1))/2;
    centreY=(axis(J(:))+axis(J(:)+1))/2;
    keep=abs(centreX)>1.05 | abs(centreY)>1.05;
    aa=I(keep)+(J(keep)-1)*N; bb=aa+1; cc=aa+N; dd=cc+1;
    elements=[aa,bb,dd;aa,dd,cc];
    used=unique(elements(:)); map=zeros(size(nodes,1),1); map(used)=1:numel(used);
    nodes=nodes(used,:); elements=reshape(map(elements),size(elements));
end
obs=fcp_mesh_data(nodes,elements);
obs.n=n;
end
