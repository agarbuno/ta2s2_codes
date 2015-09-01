%% Setting the workspace
% clc;
% clear;

rng default;
global n
n = 20; k = 2;

% func = 'limetal02non';
func = 'franke2d';
X = lhsdesign(n,k);
%ModelInfo.X = bestlh(n,k,50,20);

% Calculate observed data
for i = 1:n
    y(i,1) = feval(func, (X(i,:)));
%ModelInfo.y(i,1) = branin(ModelInfo.X(i,:));
end

meanfunc = 'h';
covfunc = 'covK';
prior = 'randn';

gamma = .5;
alpha = 0.1;
meanopt = 2;

fprintf(1, '===================================================\n');
fprintf(1, 'Problem: ................................. %s\n', func);
fprintf(1, 'Dimension ................................ %3i\n', size(X,2));
fprintf(1, 'Training runs ............................ %3i\n', n);


% rng default

%% Running aims opt
N = 2000;
lb = 10^-12;
c = 2.38/sqrt(size(X,2)+1);
fprintf(1, 'Spread parameter ......................... %1.4f\n', c);
fprintf(1, '===================================================\n\n');

% rng default

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
G = H'*C;
betahat = G\(H'*a)
m = H*betahat;

sigmahat = y'*(a - C*((H'*C)\(H'*a)))/(n-size(X,2)-2)

%% Plotting samples 
print = false;

% psiopt = min(rescale_psi);

figure(1); clf;
level_curves( min(rescale_psi) ,-5,8,30,meanopt,X, y, lb);
figure(1);
set(gca, 'FontSize', 15)
pause
hold on 
for i = 1:(size(Theta,2)/2)
    title1 = sprintf('Number of training runs = %4i', n);
    title2 = sprintf('. Temperature = %4.2f', Tvec(i));
    title(strcat(title1, title2),  'interpreter', 'latex')
    h1 = scatter(exp(Theta(:,2*i-1)), exp(Theta(:,2*i)));
    axis([10^-5 10^8 10^-5 10^8]);
    
    if print 
    % Printing options
        set(gcf, 'PaperPosition', [0.6350 6.3500 20.3200 15.2400]);
        set(gcf, 'PaperSize', [21.0000 29.7000]); 
        titlefile = sprintf('~/Dropbox/Phd/Research/Reports/Adaptive_intlik/figures/branin_wiresamp_woadap%03d.eps' ,i);
        saveas( gcf, titlefile, 'eps2c')
    end 
    pause
    if i < (size(Theta,2)/2),   delete(h1),     end;
end
scatter(exp(thetaopt(1)), exp(thetaopt(2)), 50, 'filled', 'black')
hold off

%% Plot simulator and emulator
ngrid = 80;
xz = linspace(0,1,ngrid)'; yz = xz; [Xz, Yz] = meshgrid(xz, yz);
z = [ reshape(Xz,ngrid^2,1) reshape(Yz,ngrid^2,1) ];

k = feval(covfunc, X, thetahat, 1, z);
h = feval(meanfunc, z, meanopt);

yz = h*betahat + k * (K\( y - m ) );
% yz(yz > 10^3) = 0;
Zhat = reshape(yz, size(Xz));

% True function
figure(2); clf; hold on
Z = feval(func, [reshape(Xz, ngrid^2, 1), reshape(Yz, ngrid^2, 1)]);
Z = reshape(Z, size(Xz));

% contour(Xz, Yz, Zhat, 0:10:300, '--r'); contour(Xz, Yz, Z, 0:10:300);
% contour(Xz, Yz, Zhat, [120, 500, 1000, 10000], '--r'); contour(Xz, Yz, Z, [120, 500, 1000, 10000]);
contour3(Xz, Yz, Zhat, 20, '--r'); contour3(Xz, Yz, Z, 20); 

scatter(X(:,1), X(:,2), 'filled')
title(sprintf('Number of training runs = %4i', n),  'interpreter', 'latex')
hold off
 
%% Plot the difference between simulator and emulator
 
% figure(3); clf; hold on;
% contourf(Xz, Yz, (Zhat - Z), 25);
% scatter(X(:,1), X(:,2), 'filled', 'red');
% title(sprintf('Error emulator vs simulator. Number of training runs = %4i', n),  'interpreter', 'latex')
% hold off

%% Plot the error estimates for the most probable emulator

% figure(4); clf; hold on;
% tic 
% Kchol = chol(K)'; Kk = (Kchol\k')';
% W = speye(n); W2 = speye(size(H,2));
% 
% c_star = ones(size(k,1),1) - sum(Kk*W .* conj(Kk),2);
% 
% A = (h - k * (K\H));
% Gchol = chol(G)'; GA = (Gchol\A')';
% c_star = c_star - sum(GA*W2 .* conj(GA),2);
% toc
% % c_star = ones(size(k,1),1) - diag(k * (K\k'));
% % c_star = c_star - diag(A * (G\A'));
% 
% c_star = reshape(c_star, size(Xz));
% 
% % contourf(Xz, Yz, log(sqrt(sigmahat * c_star)), 10);
% contourf(Xz, Yz, (sqrt(sigmahat * c_star)), 50);
% scatter(X(:,1), X(:,2), 'filled', 'red');
% title(sprintf('Error estimation. Number of training runs = %4i', n),  'interpreter', 'latex')
% hold off;


%% Find the modes of the density estimation
% figure(5); clf
% 
% % Percentage of canidates from the sample
% perc = 0.50;
% ksdensity(Hnew)
% [a, b, c]= ksdensity(Hnew);
% hold on 
% plot([min(Hnew),min(Hnew)],ylim, '--')
% plot(xlim, [max(a),max(a)],'--')
% hold off
% 
% limsup = max(a);
% liminf = min(a);
% treshold = (limsup + liminf)/2;
% hold on; plot(xlim, [treshold,treshold],'--'); hold off
% 
% dom = b(a >= treshold & a<= limsup);
% dominf = min(dom);
% domsup = max(dom);
% 
% tresholdnew = treshold;
% counter = 0;
% 
% while  abs(sum(Hnew <= domsup & Hnew >= dominf)/length(Hnew) - perc ) > 10^-4 && ...
%         abs(liminf - limsup) > 10^-2 || counter == 0
%     
%     if sign(sum(Hnew <= domsup & Hnew >= dominf)/length(Hnew) - perc ) == 1
%         liminf = treshold;
%         treshold = (limsup + liminf)/2;
%         dom = b(a >= treshold & a<= max(a));
%         dominf = min(dom);
%         domsup = max(dom);
%     else
%         limsup = treshold;
%         treshold = (limsup + liminf)/2;
%         dom = b(a >= treshold & a<= max(a));
%         dominf = min(dom);
%         domsup = max(dom);
%     end
%     
%     hold on; 
%     plot(xlim, [treshold,treshold],'--'); 
%     pause
%     hold off
%     counter = counter + 1;
% end
% 
% hold on; plot([domsup,domsup], ylim,'--'); 
% plot([dominf,dominf], ylim,'--'); hold off
% 
% sum(Hnew <= domsup & Hnew >= dominf)/length(Hnew)
% % 
% %% Plot the mixtures of gaussians
% % Trials
% xz = linspace(0,1,ngrid)'; yz = xz; [Xz, Yz] = meshgrid(xz, yz);
% z = [ reshape(Xz,ngrid^2,1) reshape(Yz, ngrid^2,1) ];
% h = feval(meanfunc, z, meanopt);
% Zhatmix = zeros(size(Xz));
% 
% % Drawing the best candidates from ordered values
% % [Hsort, index] = sort(Hnew);
% % thetasort = theta(index,:);
% % rescale_psisort = rescale_psi(index,:);
% % 
% % [thetasort, index] = unique(thetasort, 'rows');
% % Hsort = Hsort(index,:);
% % rescale_psisort = rescale_psisort(index,:);
% % nsamp = floor(size(Hsort,1)*.05);
% 
% 
% % Draw the best candidates from credible intervals
% [Hsort, index] = sort(Hnew);
% thetasort = theta(index,:);
% rescale_psisort = rescale_psi(index,:);
% 
% [thetasort, index] = unique(thetasort, 'rows');
% Hsort = Hsort(index,:);
% rescale_psisort = rescale_psisort(index,:);
% 
% index = Hsort <= domsup & Hsort >= dominf;
% Hsort = Hsort(index);
% thetasort = thetasort(index,:);
% rescale_psisort = rescale_psisort(index,:);
% 
% nsamp = length(Hsort)
% w_new = zeros(nsamp,1);
% 
% for i = 1:nsamp
%     w_new(i) = exp((likelihood(X, y, meanfunc, 'covK', exp(thetasort(i,:)),...
%             rescale_psisort(i,:), meanopt) - 700));
% end
% w_new = w_new./sum(w_new);
% 
% meanyz = zeros(numel(Xz), nsamp);
% varyz = zeros(numel(Xz), nsamp);
% 
% tic
% 
% for i = 1:nsamp
%     K = feval(covfunc, X, exp(thetasort(i,:)), 1);
%     K = K + rescale_psisort(i,:)  * eye(n);
% 
%     a = K\y;    C = K\H;    G = H'*C;    betahat = G\(H'*a);
%     m = H*betahat;
% 
%     sigmahat = y'*(a - C*((H'*C)\(H'*a)))/(n-size(X,2)-2); 
%     
%     k = feval(covfunc, X, exp(thetasort(i,:)), 1, z);
%     yz = h*betahat + k * (K\( y - m ) );
%     
%     Zhatmix = Zhatmix + w_new(i) * reshape(yz, size(Xz));
%     meanyz(:,i) = yz;
%     
%     Kchol = chol(K)'; Kk = (Kchol\k')';
% 
%     % This is an efficient way to compute the diagonal of a product of
%     % matrices. In comments are the past computations, more readable but
%     % way more expensive.
%     %     c_star = ones(size(k,1),1) - diag(k*(K\k'));
%     c_star = ones(size(k,1),1) - sum(Kk*W .* conj(Kk),2);
% 
%     A = (h - k * (K\H)); Gchol = chol(G)'; GA = (Gchol\A')';
%     c_star = c_star - sum(GA*W2 .* conj(GA),2);
%     %     c_star = c_star - diag(A * (G\A'));
%     
%     varyz(:,i) = sigmahat * c_star;
% end
% 
% clear K a C G A k yz Kk Gcohl W W2 GA;
% 
% toc
% 
% figure(6); clf; hold on;
% contour(Xz, Yz, Zhatmix, 0:10:300, '--r'); contour(Xz, Yz, Z, 0:10:300);
% % contour(Xz, Yz, Zhatmix, [120, 500, 1000, 10000], '--r'); contour(Xz, Yz, Z, [120, 500, 1000, 10000]);
% % contour(Xz, Yz, Zhat, '--r'); contour(Xz, Yz, Z); 
% scatter(X(:,1), X(:,2), 'filled')
% title(sprintf('Mean estimation w/mixture. Number of training runs = %4i', n),  'interpreter', 'latex')
% hold off
% 
% %% Plot the estimates of the error from the mixture of emulators
% c_starmix = zeros(size(Xz));
% tic 
% for i = 1:nsamp
%     c_starmix = c_starmix + w_new(i) * ...
%         ((reshape(meanyz(:,i), size(Xz)) - Zhatmix).^2 + ...
%         reshape(varyz(:,i), size(Xz)));
% end
% toc
% % for i = 1:nsamp
% %     c_starmix = c_starmix + w_new(i) * ...
% %         ((reshape(meanyz(:,i), size(Xz))).^2 - reshape(Zhatmix, size(Xz)).^2 + ...
% %         reshape(varyz(:,i), size(Xz)));
% % end
% 
% figure(7); clf; hold on;
% contourf(Xz, Yz, sqrt(c_starmix), 50);
% scatter(X(:,1), X(:,2), 'filled', 'red');
% title(sprintf('Error estimation w/mixture. Number of training runs = %4i', n),  'interpreter', 'latex')
% hold off;
% 
% %% Plot the estimation between the mixture and the most probable 
% 
% figure(8); clf; hold on;
% contour(Xz, Yz, Zhatmix, 0:10:300, '--r'); contour(Xz, Yz, Zhat, 0:10:300, '--b');
% contour(Xz, Yz, Zhatmix, '--r'); contour(Xz, Yz, Zhat, '--b');
% scatter(X(:,1), X(:,2), 'filled');
% title(sprintf('Mixture vs MAP. Number of training runs = %4i', n),  'interpreter', 'latex')
% hold off;
% % 
% % 
% % 
%% Plot the energy level curves. To see how does the changes in the 
% temperature affects the landscape

% L = [];
% 
% for i = 1:length(Tvec)
%     figure(5); clf;
%     L = energy_curves(psiopt,-5,8,30,meanopt,X, y, Tvec(i), L);
% %     hold on
% %         scatter(exp(Theta(:,2*i-1)), exp(Theta(:,2*i)));
% %     hold off
%     pause
% end

residuals_plots

figure(2); clf;
ylabel('Standarized Residuals', 'interpreter', 'latex')
subplot(2,1,1);
hold on
plot(res_ind, '.', 'MarkerSize',10)
title('Residual plot using MAP' ,  'interpreter', 'latex')
plot(get(gca,'xlim'), [3 3], '--k'); 
plot(get(gca,'xlim'), [-3 -3], '--k'); 
set(gca, 'FontSize', 15)
axis([0 100 -10 10])
hold off

subplot(2,1,2);
hold on
plot(res_ind_mix, '.', 'MarkerSize',10)
xlabel('Index', 'interpreter', 'latex')
title('Residual plot using mixture model' ,  'interpreter', 'latex')
plot(get(gca,'xlim'), [3 3], '--k'); 
plot(get(gca,'xlim'), [-3 -3], '--k'); 
set(gca, 'FontSize', 15)
axis([0 100 -10 10])
hold off

%% Plotting emulator and mix emulator

figure(3); clf; hold on
Z = feval(func, [reshape(Xz, ngrid^2, 1), reshape(Yz, ngrid^2, 1)]);
Z = reshape(Z, size(Xz));

contour(Xz, Yz, Z, 30)

% contour(Xz, Yz, Zhat, 0:10:300, '--r'); contour(Xz, Yz, Z, 0:10:300);
% contour(Xz, Yz, Zhat, [120, 500, 1000, 10000], '--r'); contour(Xz, Yz, Z, [120, 500, 1000, 10000]);
% contour3(Xz, Yz, Zhat, 20, '--r'); contour3(Xz, Yz, Z, 20); 

% scatter(X(:,1), X(:,2), 'filled')
% title(sprintf('Number of training runs = %4i', n),  'interpreter', 'latex')
hold off