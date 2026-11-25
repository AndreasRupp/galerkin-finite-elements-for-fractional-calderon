function x = fcp_randn(seed,n)
validateattributes(seed,{'numeric'},{'scalar','integer','nonnegative','<',1e7});
validateattributes(n,{'numeric'},{'scalar','integer','nonnegative'});
m1=2147483563; a1=40014; m2=2147483399; a2=40692; % Combined generator of L'Ecuyer, Comm. ACM 31 (1988), 742-751.
s1=mod(7*seed^2+69069*seed+1234567,m1-1)+1;
s2=mod(11*seed^2+30269*seed+7654321,m2-1)+1;
pairs=ceil(n/2); count=2*pairs; u=zeros(count,1);
for k=-49:count
    s1=mod(a1*s1,m1); s2=mod(a2*s2,m2);
    if k>0
        z=s1-s2;
        if z<1, z=z+m1-1; end
        u(k)=z/m1;
    end
end
radius=sqrt(-2*log(u(1:2:end))); angle=2*pi*u(2:2:end); % Box-Muller normal variates.
x=[radius.*cos(angle);radius.*sin(angle)];
x=x(1:n);
end
