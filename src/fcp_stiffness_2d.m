function [A,info] = fcp_stiffness_2d(mesh,order)
if nargin<2, order=8; end
if mesh.dim~=2, error('fcp:dimension','A triangular two-dimensional mesh is required.'); end
ne=size(mesh.elements,1); nn=size(mesh.nodes,1);
info.quadrature_order=order;
info.maximum_pair_difference=0;
info.unresolved_pairs=0;
info.reference_pairs=0;
info.full_space=true;
[structured,step,cells,shape,n]=structured_geometry(mesh);
if structured
    [T11,T12,pairs]=fcp_pair_table(n+1,order,true); % Triangle-pair integrals for the translation-invariant grid.
    info.maximum_pair_difference=pairs.maximum_pair_difference;
    info.unresolved_pairs=pairs.unresolved_pairs;
    info.reference_pairs=pairs.pairs;
    T={T11,T12;rot90(T12,2),T11};
    sx={[-1,0;1,0],[0,-1;0,1]}; sy={[0,0;-1,1],[-1,1;0,0]};
    t=zeros(2*n+1,2*n+1);
    for a=1:2
        for b=1:2
            kernel=conv2(rot90(sx{b},2),sx{a})+conv2(rot90(sy{b},2),sy{a});
            t=t+conv2(T{a,b},kernel,'same');
        end
    end
    t=(t+rot90(t,2))*(step/(4*pi));
    [I,J]=ind2sub([n+1,n+1],mesh.interior(:));
    position=I+(2*n+1)*J;
    A=t(bsxfun(@minus,position,position.')+(n+1)+n*(2*n+1));
    info.symmetry=norm(A-A.','fro')/max(norm(A,'fro'),realmin);
    A=(A+A.')/2;
    info.structured=true;
    if info.unresolved_pairs>0
        warning('fcp:quadrature','%d triangle pairs did not meet the requested quadrature comparison.',info.unresolved_pairs);
    end
    return
else
    H=zeros(ne,ne);
    for k=1:ne
        K=mesh.nodes(mesh.elements(k,:),:);
        for ell=1:k
            L=mesh.nodes(mesh.elements(ell,:),:);
            [entry,check]=fcp_triangle_pair(K,L,order);
            H(k,ell)=entry; H(ell,k)=entry;
            info=record_check(info,check);
        end
    end
end
row=repmat((1:ne).',1,3);
Gx=sparse(row(:),mesh.elements(:),reshape(mesh.grad(:,:,1),[],1),ne,nn);
Gy=sparse(row(:),mesh.elements(:),reshape(mesh.grad(:,:,2),[],1),ne,nn);
Gx=Gx(:,mesh.interior); Gy=Gy(:,mesh.interior);
A=(Gx.'*(H*Gx)+Gy.'*(H*Gy))/(2*pi); % Full-space fractional energy for s = 1/2.
info.symmetry=norm(A-A.','fro')/max(norm(A,'fro'),realmin);
A=full((A+A.')/2);
info.structured=structured;
if info.unresolved_pairs>0
    warning('fcp:quadrature','%d triangle pairs did not meet the requested quadrature comparison.',info.unresolved_pairs);
end
end

function info=record_check(info,check)
info.maximum_pair_difference=max(info.maximum_pair_difference,check.relative_difference);
info.unresolved_pairs=info.unresolved_pairs+double(~check.converged);
info.reference_pairs=info.reference_pairs+1;
end

function [ok,step,cells,shape,n]=structured_geometry(mesh)
ok=false; step=NaN; cells=[]; shape=[]; n=0;
if ~isfield(mesh,'n'), return; end
n=mesh.n;
origin=min(mesh.nodes,[],1); extent=max(mesh.nodes,[],1)-origin;
step=extent(1)/n;
if abs(extent(1)-extent(2))>1e-12*extent(1), return; end
ne=size(mesh.elements,1);
if ne~=2*n*n, return; end
centres=(mesh.nodes(mesh.elements(:,1),:)+mesh.nodes(mesh.elements(:,2),:)+mesh.nodes(mesh.elements(:,3),:))/3;
cells=floor(bsxfun(@minus,centres,origin)/step+1e-10);
local_centres=bsxfun(@minus,centres,origin)/step-cells;
shape=1+(local_centres(:,2)>local_centres(:,1));
canonical={sortrows([0,0;1,0;1,1]),sortrows([0,0;1,1;0,1])};
for k=1:ne
    vertices=bsxfun(@minus,mesh.nodes(mesh.elements(k,:),:),origin+step*cells(k,:))/step;
    if max(max(abs(sortrows(vertices)-canonical{shape(k)})))>1e-10, return; end
end
if any(cells(:)<0) || any(cells(:)>=n), return; end
ok=true;
end
