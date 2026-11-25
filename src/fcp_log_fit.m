function fit = fcp_log_fit(delta,error_mean)
delta=delta(:); error_mean=error_mean(:);
if numel(delta)~=numel(error_mean) || numel(delta)<3 || ...
        any(~isfinite(delta)) || any(delta<=0 | delta>=1) || ...
        any(~isfinite(error_mean)) || any(error_mean<=0)
    error('fcp:logfit','The logarithmic fit requires positive errors and noise levels between zero and one.');
end
x=log(abs(log(delta))); y=log(error_mean);
p=[ones(size(x)),x]\y;
fit.C=exp(p(1)); fit.gamma=-p(2);
fit.values=fit.C*abs(log(delta)).^(-fit.gamma);
fit.r2=1-sum((y-log(fit.values)).^2)/max(sum((y-mean(y)).^2),realmin);
fit.delta=delta;
end
