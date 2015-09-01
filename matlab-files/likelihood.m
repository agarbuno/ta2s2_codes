function L = likelihood(data, y, meanfunc, covfunc, theta, psi, meanopt, lb)

% For numerical stability it should be better to return the negative of the 
% log likelihood

[n, q] = size(data);

H = feval(meanfunc, data, meanopt);
K = feval(covfunc, data, theta, 1);
K = K + psi * speye(n);

%%%%%%%% Traditional computations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = K\y;
C = K\H;
G = H'*C;

sigma = y'*(a - C*(G\(H'*a)))/(n-q-2);
detK = det(K);

L = (n-q)/2 * log(sigma) + 0.5 * log(detK) + 0.5 * log(det(G)); 

% Adding the Jacobian to rescale the function
% z = (1-psi)/(psi-lb); 
% p = 1/(1+exp(-z));
% LJac = sum(log(theta)) + log(p) + log(1-p);
% L = L - LJac;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Cholesky computations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M = chol(K);
% a = M'\y;
% C = M'\H;
% u = C'*a;
% G = C'*C;
% 
% sigma = a'*a - u'*(G\u)/(n-q-2);
% logdetK = sum(2*log(diag(M)));
% 
% L = (n-q)/2 * log(sigma) + 0.5 * logdetK + 0.5 * log(det(G)); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Separate the components of likelihood to look at them individually
% L = (n-q)/2 * log(sigma);
% L =  0.5 * log(detK);
% L =  0.5 * log(det(H'*C));

% With other priors
% L = L - log(prior_lognormal(theta, 3));
%
% Reference prior
[I] = information_matrix(data, theta, G, C, K);
L = L - I;

% L = L - log(sqrt(det(I)));

L = L + 400;