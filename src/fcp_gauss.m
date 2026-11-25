function [x,w] = fcp_gauss(n)
validateattributes(n,{'numeric'},{'scalar','integer','positive'});
if n==1, x=0.5; w=1; return; end
k=(1:n-1).'; beta=k./sqrt(4*k.^2-1);
[V,D]=eig(diag(beta,1)+diag(beta,-1));
[x,id]=sort(diag(D)); w=(V(1,id).^2).'; x=(x+1)/2;
end
