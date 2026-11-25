function c = fcp_constant(dim,s)
if s<=0 || s>=1, error('fcp:order','The fractional order must lie in (0,1).'); end
c=4^s*s*gamma(dim/2+s)/(pi^(dim/2)*gamma(1-s));
end
