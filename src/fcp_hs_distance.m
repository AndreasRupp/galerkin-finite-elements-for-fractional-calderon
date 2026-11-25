function value = fcp_hs_distance(hs,v)
d=hs.Pc*v(:)-hs.ref;
value=d'*(hs.S*d);
if isfield(hs,'offset2'), value=value+hs.offset2; end
value=sqrt(max(0,value));
end
