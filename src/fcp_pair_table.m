function [T11,T12,info] = fcp_pair_table(n,order,compare)
if nargin<2, order=6; end
if nargin<3, compare=false; end
validateattributes(n,{'numeric'},{'scalar','integer','>=',1});
K1=[0,0;1,0;1,1]; K2=[0,0;1,1;0,1];
width=2*n-1; offsets=-(n-1):(n-1);
[DX,DY]=ndgrid(offsets,offsets);
close_pairs=max(abs(DX),abs(DY))<=5;
info.maximum_pair_difference=0; info.unresolved_pairs=0; info.pairs=0;
info.far_order=order+2; info.compared=compare;

T11=zeros(width,width);
half=DX>0 | (DX==0 & DY>=0);
ids=find(half & close_pairs);
for k=1:numel(ids)
    [T11(ids(k)),check]=fcp_triangle_pair(bsxfun(@plus,K1,[DX(ids(k)),DY(ids(k))]),K1,order);
    info=record(info,check);
end
ids=find(half & ~close_pairs);
[T11(ids),info]=far_pairs([DX(ids),DY(ids)],K1,K1,order,compare,info);
mirror=rot90(T11,2);
T11(~half)=mirror(~half);

T12=zeros(width,width);
ids=find(close_pairs);
for k=1:numel(ids)
    [T12(ids(k)),check]=fcp_triangle_pair(bsxfun(@plus,K1,[DX(ids(k)),DY(ids(k))]),K2,order);
    info=record(info,check);
end
ids=find(~close_pairs);
[T12(ids),info]=far_pairs([DX(ids),DY(ids)],K1,K2,order,compare,info);
end

function [values,info]=far_pairs(shift,K,L,order,compare,info)
values=far_values(shift,K,L,order+2);
info.pairs=info.pairs+numel(values);
if compare && ~isempty(values)
    lower=far_values(shift,K,L,order);
    info.maximum_pair_difference=max(info.maximum_pair_difference,max(abs(values-lower)./abs(values)));
end
if any(~isfinite(values) | values<=0)
    error('fcp:triangle','A positive triangle-pair integral was not resolved.');
end
end

function info=record(info,check)
info.maximum_pair_difference=max(info.maximum_pair_difference,check.relative_difference);
info.unresolved_pairs=info.unresolved_pairs+double(~check.converged);
info.pairs=info.pairs+1;
end

function values=far_values(shift,K,L,q)
[xk,wk]=triangle_rule(K,q); [xl,wl]=triangle_rule(L,q);
count=size(shift,1); values=zeros(count,1); np=numel(wk);
block=max(1,floor(1.5e6/np));
for first=1:block:count
    ids=(first:min(first+block-1,count)).'; m=numel(ids);
    px=bsxfun(@plus,shift(ids,1),xk(:,1).'); py=bsxfun(@plus,shift(ids,2),xk(:,2).');
    first_integral=reshape(fcp_triangle_potential([px(:),py(:)],L),m,np)*wk;
    px=bsxfun(@minus,xl(:,1).',shift(ids,1)); py=bsxfun(@minus,xl(:,2).',shift(ids,2));
    second_integral=reshape(fcp_triangle_potential([px(:),py(:)],K),m,np)*wl;
    values(ids)=(first_integral+second_integral)/2;
end
end

function [points,weights]=triangle_rule(K,q)
[t,w]=fcp_gauss(q);
[U,V]=ndgrid(t,t); weights=(w*w.').*(1-U);
points=bsxfun(@plus,K(1,:),U(:)*(K(2,:)-K(1,:))+((1-U(:)).*V(:))*(K(3,:)-K(1,:)));
weights=abs(det([K(2,:)-K(1,:);K(3,:)-K(1,:)]))*weights(:);
end
