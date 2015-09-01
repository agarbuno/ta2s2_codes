function z=bnnLikelihood(theta,x,y,n)

sigma_prior = 5;

for i=1:n
    mu(i)=y(i)-theta(1)*activation(x(i)*theta(3)+ theta(5))-theta(2)*activation(x(i)*theta(4)+theta(6));
    a(i)=normpdf(mu(i),0,exp(-theta(7)/2));
end
z=-sum(log(a));
z = z + prior(theta,sigma_prior) + 100;
