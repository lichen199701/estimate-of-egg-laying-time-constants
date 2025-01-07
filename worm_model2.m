%   Probability density function for intervals 
%   3-state model
%   This function evaluates Eq. (12) of the paper, by Zhou, Schafer and 
%   Schafer, IEEE Trans. Signal Processing, vol. 46, no. 10, October 1998.  

function f=worm_model2(x,lambda3,lambda1,lambda2)

p=lambda3/(lambda1+lambda3);
k1=p*(lambda1-lambda2)/(lambda1-p*lambda2);     % weight
k2=lambda1*(1-p)/(lambda1-p*lambda2);           % weight

f=k1*lambda1*exp(-lambda1*x)+k2*p*lambda2*exp(-p*lambda2*x);       
end

