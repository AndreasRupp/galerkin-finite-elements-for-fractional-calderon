function P = fcp_evaluate(mesh,points)
if mesh.dim==1, points=points(:); end
np=size(points,1); nv=mesh.dim+1;
rows=[]; cols=[]; vals=[];
if mesh.dim==1
    x=mesh.nodes(:,1); tolerance=64*eps*max(1,max(abs(x)));
    active=find(points(:,1)>=x(1)-tolerance & points(:,1)<=x(end)+tolerance);
    if isempty(active), P=sparse(np,numel(mesh.interior)); return; end
    xp=min(max(points(active,1),x(1)),x(end));
    [~,cell]=histc(xp,x); cell=min(cell,numel(x)-1);
    local=(xp-x(cell))./(x(cell+1)-x(cell));
    nodes=[cell,cell+1]; lambda=[1-local,local];
else
    x=unique(mesh.nodes(:,1)); y=unique(mesh.nodes(:,2));
    nx=numel(x)-1; ny=numel(y)-1;
    if size(mesh.nodes,1)~=(nx+1)*(ny+1)
        error('fcp:evaluate','Two-dimensional evaluation requires the full rectangular mesh.');
    end
    hx=(x(end)-x(1))/nx; hy=(y(end)-y(1))/ny;
    if max(abs(diff(x)-hx))>128*eps || max(abs(diff(y)-hy))>128*eps
        error('fcp:evaluate','Two-dimensional evaluation requires a uniform grid.');
    end
    tolerance=64*eps*max(1,max(abs(mesh.nodes(:))));
    active=find(points(:,1)>=x(1)-tolerance & points(:,1)<=x(end)+tolerance & points(:,2)>=y(1)-tolerance & points(:,2)<=y(end)+tolerance);
    if isempty(active), P=sparse(np,numel(mesh.interior)); return; end
    xp=min(max(points(active,1),x(1)),x(end)); yp=min(max(points(active,2),y(1)),y(end));
    ix=min(floor((xp-x(1))/hx),nx-1); iy=min(floor((yp-y(1))/hy),ny-1);
    rx=(xp-(x(1)+ix*hx))/hx; ry=(yp-(y(1)+iy*hy))/hy;
    a=ix+1+iy*(nx+1); b=a+1; c=a+(nx+1); d=c+1;
    lower=ry<=rx; nodes=[a,b,d]; lambda=[1-rx,rx-ry,ry];
    nodes(~lower,:)=[a(~lower),d(~lower),c(~lower)];
    lambda(~lower,:)=[1-ry(~lower),rx(~lower),ry(~lower)-rx(~lower)];
end
for j=1:nv
    col=mesh.full_to_interior(nodes(:,j)); take=col>0;
    rows=[rows;active(take)]; cols=[cols;col(take)]; vals=[vals;lambda(take,j)];
end
P=sparse(rows,cols,vals,np,numel(mesh.interior));
end
