function b = fcp_exterior_load_exact(mesh,datum,s,order)
[y,wy]=fcp_datum_quadrature(datum,mesh.dim);
weights=wy.*fcp_exterior_datum(y,datum);
keep=weights~=0; y=y(keep,:); weights=weights(keep);
if mesh.dim==2
    if abs(s-0.5)>64*eps
        error('fcp:order','The analytic two-dimensional load is implemented for s=1/2.');
    end
    b=zeros(numel(mesh.interior),1);
    for first=1:256:size(y,1)
        ids=first:min(first+255,size(y,1));
        b=b+fcp_observation_values_2d(mesh,y(ids,:))'*weights(ids);
    end
    return
end
[x,wx,Phi]=fcp_quadrature(mesh,order);
load=zeros(size(wx)); c=fcp_constant(1,s);
for first=1:512:size(x,1)
    ids=first:min(first+511,size(x,1));
    r=abs(bsxfun(@minus,x(ids),y'));
    if any(r(:)<=0), error('fcp:separation','The exterior datum must be separated from Omega.'); end
    load(ids)=-c*(r.^(-(1+2*s)))*weights;
end
b=Phi'*(wx.*load);
end
