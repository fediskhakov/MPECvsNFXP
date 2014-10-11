function [P, y]=discretized_normal(y,  mu, sigma)
	n=numel(y);
	d=(y(n)-y(1))/(n-1);
	p_k= @ (k) normcdf((y(k) + d-mu)/sigma) - normcdf((y(k) -mu)/sigma);
	P=nan(n,1);
	P(1)= normcdf((y(1) + d-mu)/sigma);
	P(n)= 1 - normcdf((y(n) -mu)/sigma);
	for k=2:n-1
		P(k)=p_k(k);
	end
end