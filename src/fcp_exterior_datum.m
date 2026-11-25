function f = fcp_exterior_datum(x,datum)
f=zeros(size(x,1),1);
for j=1:size(datum.centres,1)
    r2=sum(bsxfun(@minus,x,datum.centres(j,:)).^2,2)/datum.radius^2;
    active=r2<1;
    f(active)=f(active)+exp(datum.offset-1./(1-r2(active)));
end
end
