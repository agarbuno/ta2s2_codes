%% Setting the workspace
% clc;
clear;

rng default;
global n
n = 500;

func = 'wingweight';
[X,y] = startup(n, func);
[Xtry,Ytry] = startup(500, func);

meanfunc = 'h';
covfunc = 'covK';
prior = 'randn';

gamma = .3;
alpha = 0.05;
meanopt = 1;

%% Running aims opt
N = 5000;
lb = 10^-12;
c = 2.38/sqrt(size(X,2)+1);

fprintf(1, '===================================================\n');
fprintf(1, 'Problem: ................................. %s\n', 'Wing Weight');
fprintf(1, 'Dimension ................................ %3i\n', size(X,2));
fprintf(1, 'Training runs ............................ %3i\n', size(X,1));
fprintf(1, 'Spread parameter ......................... %1.4f\n', c);
fprintf(1, '===================================================\n\n');

% rng default;

tic
    [theta, psinew, Hnew, k, w, Theta, Accep, Tvec] = slice_parallel_aims_opt(X, y, gamma, alpha, N, c, ...
        meanfunc, covfunc, meanopt, prior, lb);
toc

fprintf(1, '\nInterval: [ %4.8f, %4.8f ] \n', min(Hnew), max(Hnew));


%% Nugget term 
rescale_psi = (1-lb)./(1+exp(-psinew)) + lb;
psiopt = rescale_psi(Hnew == min(Hnew));
psiopt = unique(psiopt);

thetaopt = theta(Hnew == min(Hnew),:);
thetaopt = unique(thetaopt,'rows');

fprintf(1, 'Optimal nugget ............... %8.2e\n', psiopt);
if size(thetaopt,2) == 2
    fprintf(1, 'Optimal theta: [ %4.8f, %4.8f ] \n', exp(thetaopt(1)), exp(thetaopt(2)) );
else
    fprintf(1, 'Optimal theta: [ %4.8f, %4.8f, %4.8f, %4.8f, %4.8f ] \n', ...
        exp(thetaopt(1)/2), exp(thetaopt(2)/2), exp(thetaopt(3)/2), exp(thetaopt(4)/2), exp(thetaopt(5)/2) );
end

exp(-thetaopt)

[~, index ] = sort(-exp(-thetaopt))

% Full regression model
%      7     8     4     3     1     9     6     5     2    10
%      8     7     3     4     1     9     6     5    10     2
%      7     8     3     4     1     9     6    10     5     2

% Constant mean model
%      8     7     3     1     9     4     6    10     5     2


%% Computes estimations for sigma
params = [exp(theta)];
L = -Hnew;
n = size(X,1);

maxp = params(L == max(L),:);
Lmax = max(L);
maxp = maxp(1,:);

thetahat = maxp;

H = feval(meanfunc, X, meanopt);
K = feval(covfunc, X, thetahat, 1);
K = K + psiopt * eye(n);

a = K\y;
C = K\H;
betahat = (H'*C)\(H'*a)
m = H*betahat;

sigmahat = y'*(a - C*((H'*C)\(H'*a)))/(n-size(X,2)-2)

%% Computes estimations for sigma
params = exp(theta);
L = -Hnew;
n = size(X,1);

maxp = params(L == max(L),:);
Lmax = max(L);
maxp = maxp(1,:);

thetahat = maxp;

H = feval(meanfunc, X, meanopt);
K = feval(covfunc, X, thetahat, 1);
K = K + psiopt * eye(n);

a = K\y;    C = K\H;    G = H'*C;
betahat = G\(H'*a);
m = H*betahat;

sigmahat = y'*(a - C*((H'*C)\(H'*a)))/(n-size(X,2)-2);

%% Plot standarized residuals

k = feval(covfunc, X, thetahat, 1, Xtry);
h = feval(meanfunc, Xtry, meanopt);

% Estimated mean
mean_y = h*betahat + k * (K\( y - m ) );
RMSE = sqrt(mean((mean_y - Ytry).^2))


% Error prediction
Kchol = chol(K)'; Kk = (Kchol\k')';
c_star = ones(size(k,1),1) - sum(Kk .* conj(Kk),2);
A = (h - k * (K\H)); Gchol = chol(G)'; GA = (Gchol\A')';
c_star = c_star - sum(GA .* conj(GA),2);
var_y = sigmahat * c_star;

% Standarized independent residuals
res_ind = (Ytry - mean_y)./sqrt(var_y);

% Plotting the residuals
% Plotting the residuals
figure(1); clf;
ylabel('Standarized Residuals', 'interpreter', 'latex')
subplot(2,1,1);
hold on
plot(res_ind, '.', 'MarkerSize',10)
title('Residual plot using MAP' ,  'interpreter', 'latex')
plot(get(gca,'xlim'), [3 3], '--k'); 
plot(get(gca,'xlim'), [-3 -3], '--k'); 
set(gca, 'FontSize', 15)
axis([0 Inf -10 10])
hold off

%plot(mean_y, res_ind, '.')
%plot(Xtry(:,1), res_ind,'.')
%plot(Xtry(:,2), res_ind,'.')

%% Using the mixture

y_mix = zeros(size(Xtry,1),1);
w_new = 1/size(theta,1) * ones(size(theta,1),1);
meanyz = zeros(size(Xtry,1), size(theta,1));
varyz = zeros(size(Xtry,1), size(theta,1));

% This for iterates through all the samples from paims
for i = 1:size(theta,1)
    K = feval(covfunc, X, exp(theta(i,:)), 1);
    K = K + rescale_psi(i,:) * speye(n);

    a = K\y;    C = K\H;    G = H'*C;    betahat = G\(H'*a);
    m = H*betahat;

    sigmahat = y'*(a - C*((H'*C)\(H'*a)))/(n-size(X,2)-2); 
    
    k = feval(covfunc, X, exp(theta(i,:)), 1, Xtry);
    yz = h*betahat + k * (K\( y - m ) );
    
    y_mix = y_mix + w_new(i) * yz;
    meanyz(:,i) = yz;
    
    Kchol = chol(K)'; Kk = (Kchol\k')';

    % This is an efficient way to compute the diagonal of a product of
    % matrices. In comments are the past computations, more readable but
    % way more expensive.
    %     c_star = ones(size(k,1),1) - diag(k*(K\k'));
    c_star = ones(size(k,1),1) - sum(Kk .* conj(Kk),2);

    A = (h - k * (K\H)); Gchol = chol(G)'; GA = (Gchol\A')';
    c_star = c_star - sum(GA .* conj(GA),2);
       
    varyz(:,i) = sigmahat * c_star;
end

RMSE = sqrt(mean((y_mix - Ytry).^2))
var_mix = bsxfun(@plus, bsxfun(@minus, meanyz, y_mix).^2, varyz) * w_new;

% Standarized independent residuals
res_ind_mix = (Ytry - y_mix)./sqrt(var_mix);

% Plotting the residuals
subplot(2,1,2);
hold on
plot(res_ind_mix, '.', 'MarkerSize',10)
xlabel('Index', 'interpreter', 'latex')
title('Residual plot using mixture model' ,  'interpreter', 'latex')
plot(get(gca,'xlim'), [3 3], '--k'); 
plot(get(gca,'xlim'), [-3 -3], '--k'); 
set(gca, 'FontSize', 15)
axis([0 Inf -10 10])
hold off

figure(2); clf; 
boxplot(theta)
set(gca, 'FontSize', 15)
