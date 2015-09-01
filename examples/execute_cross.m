%% Setting the workspace
% clc;
clear;

% rng default;
global n
dim = 2;
a = 10;

gamma = .5;
alpha = 0.01;
meanopt = 2;

% rng default

%% Running aims opt
func = 'objcross';
fprintf(1, '===================================================\n');
fprintf(1, 'Problem: ................................. %s\n', func);
fprintf(1, 'Dimension ................................ %3i\n', dim);
N = 2000;
c = 2.38/sqrt(dim+1);
fprintf(1, 'Spread parameter ......................... %1.4f\n', c);
fprintf(1, '===================================================\n\n');

% rng default

tic
    [theta, Hnew, k, w, Theta, Accep, Tvec] = slice_opt(func, dim, gamma, alpha, N, c, a);
toc

fprintf(1, '\nInterval: [ %4.8f, %4.8f ] \n', min(Hnew), max(Hnew));

figure(1); clf; 
for i = 1:(size(Theta,2)/2)
    transparentScatter(Theta(:,2*i-1), Theta(:,2*i), 0.1, (i+eps)/(k+1));
    axis([0 10 0 10])
end

transparentScatter(theta(:,1), theta(:,2), 0.1, 1, 1);

%% Running aims opt
func = 'objcenter';
fprintf(1, '===================================================\n');
fprintf(1, 'Problem: ................................. %s\n', func);
fprintf(1, 'Dimension ................................ %3i\n', dim);
c = 2.38/sqrt(dim+1);
fprintf(1, 'Spread parameter ......................... %1.4f\n', c);
fprintf(1, '===================================================\n\n');

% rng default

tic
    [theta, Hnew, k, w, Theta, Accep, Tvec] = slice_opt(func, dim, gamma, alpha, N, c, a);
toc

fprintf(1, '\nInterval: [ %4.8f, %4.8f ] \n', min(Hnew), max(Hnew));

figure(2); clf; 
for i = 1:(size(Theta,2)/2)
    transparentScatter(Theta(:,2*i-1), Theta(:,2*i), 0.1, (i+eps)/(k+1));
    axis([0 10 0 10])
end

transparentScatter(theta(:,1), theta(:,2), 0.1, 1, 1);

%% Running aims opt
alpha = 0.05;
func = 'objcorners';
fprintf(1, '===================================================\n');
fprintf(1, 'Problem: ................................. %s\n', func);
fprintf(1, 'Dimension ................................ %3i\n', dim);
c = 2.38/sqrt(dim+1);
fprintf(1, 'Spread parameter ......................... %1.4f\n', c);
fprintf(1, '===================================================\n\n');

% rng default

tic
    [theta, Hnew, k, w, Theta, Accep, Tvec] = slice_opt(func, dim, gamma, alpha, N, c, a);
toc

fprintf(1, '\nInterval: [ %4.8f, %4.8f ] \n', min(Hnew), max(Hnew));

figure(3); clf; 
for i = 1:(size(Theta,2)/2)
    transparentScatter(Theta(:,2*i-1), Theta(:,2*i), 0.1, (i+eps)/(k+1));
    axis([0 10 0 10])
end

transparentScatter(theta(:,1), theta(:,2), 0.1, 1, 1);


