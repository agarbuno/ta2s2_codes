function z=prior(theta,sigma)


z1=-log(normpdf(theta(1),0,sigma));
z2=-log(normpdf(theta(2),0,sigma));
z3=-log(normpdf(theta(3),0,sigma));
z4=-log(normpdf(theta(4),0,sigma));
z5=-log(normpdf(theta(5),0,sigma));
z6=-log(normpdf(theta(6),0,sigma));
z7=-log(normpdf(theta(7),0,sigma));

z=z1+z2+z3+z4+z5+z6+z7;
