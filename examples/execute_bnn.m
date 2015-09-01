%% Set working directory and parameters
% Add path with source files for the bayes neural network example
addpath('./bayesnn/')

% Prepare the synthetic data set
% clc; 
rng default
BNN

% Parameters for algorithm
N = 2000;
gamma = .3;
alpha = 0.01;
dim = size(thetaOPT,2);
c = 2.38/sqrt(dim+1);

%% Running paass
func = 'Bayes Neural Network';
fprintf(1, '=================================================================\n');
fprintf(1, 'Problem: ................................. %s\n', func);
fprintf(1, 'Dimension ................................ %3i\n', size(thetaOPT,2));
fprintf(1, 'Training set ............................. %3i\n', size(x,2));
fprintf(1, 'Spread parameter ......................... %1.4f\n', c);
fprintf(1, 'Samples produced ......................... %4i\n', N);
fprintf(1, '=================================================================\n\n');

tic
    [theta, Hnew, k, w, Theta, Accep, Tvec] = slice_opt_bnn(x, y, dim, gamma, alpha, N, c);
toc

%% Presenting the results
fprintf(1, '\nInterval: [ %4.8f, %4.8f ] \n', min(Hnew), max(Hnew));
fprintf(1, 'Optimal theta: '); 
theta((Hnew == min(Hnew)),:)

% Cycling through the theta samples
ypred = zeros(size(theta,1), length(x));
for i = 1:size(theta,1)
    ypred(i,:) = pred_neural_net(theta(i,:), x); 
end

figure(1); clf;
plot(x,y, 'LineWidth', 2); hold on
plot(x, mean(ypred), 'LineWidth', 2);
plot(x, mean(ypred) + 1.96 * std(ypred), '--k')
plot(x, mean(ypred) - 1.96 * std(ypred), '--k')
title('Fitted vs actual model', 'interpreter', 'latex')
legend('True', 'Fitted')
set(gca, 'FontSize', 15)
hold off

figure(2); clf;
hist(theta(:,1));
title('Estimated density $\theta_1$', 'interpreter', 'latex')
set(gca, 'FontSize', 15)

figure(3); clf;
hist(theta(:,2));
title('Estimated density $\theta_2$', 'interpreter', 'latex')
set(gca, 'FontSize', 15)

figure(4); clf;
hist(theta(:,3));
title('Estimated density $\theta_3$', 'interpreter', 'latex')
set(gca, 'FontSize', 15)

figure(5); clf;
hist(theta(:,4));
title('Estimated density $\theta_4$', 'interpreter', 'latex')
set(gca, 'FontSize', 15)

figure(6); clf;
hist(theta(:,5));
title('Estimated density $\theta_5$', 'interpreter', 'latex')
set(gca, 'FontSize', 15)

figure(7); clf;
hist(theta(:,6));
title('Estimated density $\theta_6$', 'interpreter', 'latex')
set(gca, 'FontSize', 15)

figure(8); clf;
hist(exp(-theta(:,7)/2));
title('Estimated density $\theta_7$', 'interpreter', 'latex')
set(gca, 'FontSize', 15)