%  Example: Bayesian Neural Network (BNN)
clear;
%------------------------- DATA: (x,y) ------------------------------
n=100;     % number of data points
alpha1=5;  % true values of parameters 
alpha2=-5;
beta1=-1;
beta2=-3;
gamma1=5;
gamma2=2;
sigma_eps=0.1; % std of the pediction-error
epsilon=randn(1,n)*sigma_eps;  % prediction-error
x=0.1:0.1:10;
for i=1:n
    y(i)=alpha1*activation(beta1*x(i)+gamma1)+alpha2*activation(beta2*x(i)+gamma2)+epsilon(i);
end
sigma_prior=5; % std of the prior
%----------------------------------------------------------------------
thetaOPT=[alpha1,alpha2,beta1,beta2,gamma1,gamma2,log(sigma_eps^(-2))]
% theta=randn(1,7)*sigma_prior;       % a point in the parameter space
Ptheta=prior(thetaOPT,sigma_prior)    % value of the prior pdf at theta
Ltheta=bnnLikelihood(thetaOPT,x,y,n)     % value of the likelihood at theta