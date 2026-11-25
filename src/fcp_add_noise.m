function [noisy,absolute,noise] = fcp_add_noise(model,mu,relative,seed,generator)
if nargin<5, generator='portable'; end
validateattributes(relative,{'numeric'},{'scalar','real','finite','nonnegative'});
mu=mu(:);
if strcmp(generator,'portable')
    xi=fcp_randn(seed,numel(mu));
elseif strcmp(generator,'native')
    rng(seed,'twister'); xi=randn(size(mu));
else
    error('fcp:noise','The noise generator must be portable or native.');
end
mu_norm=sqrt(max(0,mu'*(model.W*mu)));
xi_norm=sqrt(max(0,xi'*(model.W*xi)));
absolute=relative*mu_norm;
if absolute==0
    noise=zeros(size(mu));
else
    noise=(absolute/xi_norm)*xi;
end
noisy=mu+noise;
end
