function [vertices,faces,values] = fcp_surface_jumps(mesh,c,n)
if nargin<3, n=1; end
validateattributes(n,{'numeric'},{'scalar','integer','positive'});
ne=size(mesh.elements,1);
if mesh.dim~=2 || ~isequal(size(c),[ne,3]) || any(~isfinite(c(:)))
    error('fcp:surfaceJumps','The surface requires finite values at the three vertices of each triangle.');
end
local=[1,2;2,3;3,1];
ends=zeros(3*ne,2); traces=zeros(3*ne,2);
for k=1:3
    rows=(k-1)*ne+(1:ne);
    ends(rows,:)=mesh.elements(:,local(k,:));
    traces(rows,:)=c(:,local(k,:));
end
reversed=ends(:,1)>ends(:,2);
ends(reversed,:)=ends(reversed,[2,1]);
traces(reversed,:)=traces(reversed,[2,1]);
[ends,order]=sortrows(ends,[1,2]); traces=traces(order,:);
first=[1;find(any(diff(ends,1,1),2))+1]; last=[first(2:end)-1;size(ends,1)];
vertices=zeros(0,2); values=zeros(0,1); faces=zeros(0,3);
for k=1:numel(first)
    a=first(k); b=last(k);
    if b-a>1, error('fcp:surfaceMesh','An edge has more than two adjacent triangles.'); end
    ids=ends(a,:); xy=mesh.nodes(ids,:);
    t=linspace(0,1,n+1).';
    points=(1-t)*xy(1,:)+t*xy(2,:);
    left=(1-t)*traces(a,1)+t*traces(a,2);
    if a==b, right=zeros(size(t));
    else, right=(1-t)*traces(b,1)+t*traces(b,2); end
    for j=1:n
        xyj=points(j:j+1,:); lj=left(j:j+1); rj=right(j:j+1);
        difference=lj-rj;
        cuts=[0;1];
        if difference(1)*difference(2)<0
            cuts=[0;difference(1)/(difference(1)-difference(2));1]; % Split crossing traces before drawing the jump.
        end
        for m=1:numel(cuts)-1
            tau=cuts(m:m+1);
            p=(1-tau)*xyj(1,:)+tau*xyj(2,:);
            zl=(1-tau)*lj(1)+tau*lj(2);
            zr=(1-tau)*rj(1)+tau*rj(2);
            if all(zl==zr), continue; end
            offset=size(vertices,1);
            vertices=[vertices;p(1,:);p(2,:);p(2,:);p(1,:)];
            values=[values;zl(1);zl(2);zr(2);zr(1)];
            faces=[faces;offset+[1,2,3];offset+[1,3,4]];
        end
    end
end
end
