% Program to estimate the parameters of the 3-state model for egg laying
% in C. elegans.  Model as described in Zhou, Schafer, and Schafer, IEEE
% Trans. on Signal Processing, vol. 46, No. 10, October 1998.
%
% Given the observed intervals, estimate model parameters using
% the ML algorithm. Initial estimates are provided by matching
% the (bimodal) histogram of the log intervals. 
%
% INPUT:  x: intervals between egg laying events
% OUTPUT: lambda3: rate constant of A->I
%   lambda1: rate constant of A->E
%   lambda2: time constant of I->A
%
% Usage: [lambda3,lambda1,lambda2]=ml_est_only(x) 
% SEE ALSO: local_mm.m, maxn.m, fy_fun2.m, fmins.m

%           Original program by G. Tong Zhou, Georgia Tech
%           May 20, 1997
%           Revised 5/6/99
%           Revised by RWS, 12/03/2016
%
function [lambda3,lambda1,lambda2]=ml_est_only(x) 
N=length(x); % data length
delta=.5;   %change this to modify histogram if not bimodal
[ny,y]=hist(log(x),[0:delta:10]);           % histogram of log(x)
fy=ny/(delta*N);                     % convert histogram to probabilities
yi=[min(y):max(y)/500:max(y)];       % 501 points interpolation
fyi=interp1(y,fy,yi,'PCHIP');        % fyi vs. yi smoother estimate
 
% ******************* Step 1: Initial Estimates **********************
%
[maxy,max_loc]=local_mm(yi,fyi);        % find local maxima of fyi
[foo,index2]=maxn(maxy,2);              % find two greatest local maxima
 
max_loc=max_loc(index2);                % -log(lambda1) and -log(p*lambda2)
maxy=maxy(index2);                      % k1/e and k2/e

[foo,index]=sort(max_loc);              % sort location in ascending order
max_loc=max_loc(index);                 % location in ascending order
maxy=maxy(index);                       % corresponding peak magnitude

k1_hat=exp(1)*maxy(1);                  % estimate for k1
k2_hat=exp(1)*maxy(2);                  % estimate for k2
 
lambda1_hat=1/exp(max_loc(1));          % estimate for lambda1
plambda2_hat=1/exp(max_loc(2));         % estimate for p*lambda2
 
plambda1_hat=plambda2_hat+k1_hat*(lambda1_hat-plambda2_hat);    % Eq. (36)
p_hat=plambda1_hat/lambda1_hat; lambda2_hat=plambda2_hat/p_hat; % Eq. (37)
lambda3_hat=lambda1_hat*p_hat/(1-p_hat);
para1=[lambda3_hat,lambda1_hat,lambda2_hat];  % starting values for NLLS est.
%
% ******** Step 2: Nonlinear Least-Squares Estimates *****************
%
options(14)=100*5;                      % set max number of iterations
para2=fmins('fy_fun2',para1,options,[],y,fy);   %NLLS estimate of parameters
% para2 is starting estimate for ML estimate of parameters.
%
% ******************* Step 3: ML Estimates ***********************
%
worm_pdf=@(x,lambda3,lambda1,lambda2) worm_model2(x,lambda3,lambda1,lambda2);
options = statset('MaxIter',400, 'MaxFunEvals',800);
[para3,CI3]=mle(x,'pdf',worm_pdf,'start',para2,...
    'lower',[0,0,0],'upper',[1,Inf,Inf],'options',options);  
lambda3=para3(1); lambda1=para3(2); lambda2=para3(3);  % ML estimate of parameters
  