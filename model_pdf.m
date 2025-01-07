function z=model_pdf(para)

global INTERVALS

p=para(1);
lambda1=para(2);
lambda2=para(3);

k1=p*(lambda1-lambda2)/(lambda1-p*lambda2);     % weight
k2=lambda1*(1-p)/(lambda1-p*lambda2);           % weight

z=-prod(10*k1*lambda1*exp(-lambda1*INTERVALS) ...
       +10*k2*p*lambda2*exp(-p*lambda2*INTERVALS));        
        % -likelihood function: objective function to minimize
