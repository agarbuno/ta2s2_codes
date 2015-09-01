function y = pred_neural_net(theta, x)
    n = length(x);
    epsilon=randn(1,n)*exp(-theta(7)/2); 
    y = zeros(1,n);
    for i=1:n
        y(i)=theta(1)*activation(theta(3)*x(i)+theta(5))+theta(2)*activation(theta(4)*x(i)+theta(6))+epsilon(i);
    end
end