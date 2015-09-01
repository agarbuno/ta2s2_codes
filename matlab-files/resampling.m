function [theta, psi, H, w, ratenew] = resampling(H, Tempold, Tempnew, N, gamma, theta, psi, ...
    X, y, meanfunc, covfunc, c0, rate, level, lb, meanopt, accep)

w = exp ( - H * (1/Tempnew - 1/Tempold) );
w = w / sum(w);

% Resampling by weights.
% if 1/sum(w.^2) < gamma * N
%     fprintf(1, 'Resampling needed\n');
%     indexes = randsample(length(w), N, true, w);
%     theta = theta(indexes,:); psi = psi(indexes, :); 
%     H = H(indexes);
%     w = ones(length(H),1) ./ length(H);
% end

% Resampling with gaussian noise.
% if 1/sum(w.^2) < gamma * N

% Adapting the spread
ratenew = rate;
c = c0 * ratenew.^(level);
% Checkpoints to asses good mixing due to the spread parameter
checkpoints = 50;
batches = 0;
batch_accep = 0;
n_batch = 0;

resamp_accep = 0;
if accep < 0.20
    fprintf(1, 'Resampling needed *********************************\n');

    % Preallocating memory
    thetanew = zeros(size(theta));
    psinew = zeros(size(psi));
    Hnew = zeros(N,1);
%     ratenew = zeros(size(rate));
    c = c0 * ratenew.^(1*level);
    Ntheta = size(theta,2);
    
    den = exp( -H / Tempnew );
    
    params = [theta, psi];
    % Compute the covariance matrix
%     C = weightedcov(params, w);
%     % C = weightedcov(params, (1./H)/(sum(1./H)));
%     try G = chol(C)';
%     catch
%         fprintf(1, '... Resampling - Actually just use the variances ...\n');
%         G = diag(std(params));
%     end
    G = diag(std(params));
    % G = 5 * eye(size(params,2));
    
    % Move the samples with gaussian moves.
    parfor i = 1:N
        % This computes the update of the rate of reduction for the spread
        % parameter in the adaptive steps
%         if mod(i, checkpoints) == 0
%             % batch counter
%             batches = batches + 1;
%             % delta factor 
%             delta = min(0.1, i^-.5);
%             % check whether to increase or decrease rate
%             if (resamp_accep-batch_accep)/(N - n_batch) < .20
%                 ratenew = exp(log(rate) - (Tempnew - 1) * delta);
%             else if (resamp_accep-batch_accep)/(i - n_batch) > .23
%                     ratenew = exp(log(rate) + (Tempnew - 1) * delta);
%                 end
%             end
%             c = c0 * ratenew.^(level);
%             batch_accep = resamp_accep; n_batch = i;
%         end
        
        % Center the movement in the previous particle
        m = [theta(i,:), psi(i, 1)];
        newparams = m + c * ( G * randn(Ntheta + 1, 1))';
%         ratenew(i) = rate(i);

        % Extraction of the parameters
        thetanew(i,:) = newparams(1:Ntheta);
          psinew(i,:) = newparams(Ntheta+1);
        rescale_psi = (1-lb)/(1+exp(-psinew(i,:))) + lb;
        
        
        % Recalculation of the objective function
        if Tempnew < 0
            Hnew(i) = likelihood_eig(X, y, meanfunc, covfunc, exp(thetanew(i,:)),...
                rescale_psi, meanopt);
        else
            Hnew(i) = likelihood(X, y, meanfunc, covfunc, exp(thetanew(i,:)),...
                rescale_psi, meanopt);
        end
        
        % Extraction of 
        c1 = min(1,exp(-Hnew(i)/Tempnew)/ den(i));
        
        if  (Hnew(i) == Inf || Hnew(i) == -Inf || imag(Hnew(i)) ~= 0 )  || rand(1,1) >= c1
            thetanew(i,:) = theta(i,:);
            psinew(i,:) = psi(i,:);
%             ratenew(i,1) = rate(i);
            
            Hnew(i) = H(i,:);
        else
            resamp_accep = resamp_accep + 1;
        end
        
    end
    
    theta = thetanew; psi = psinew; H = Hnew;
    w = exp ( - H * (1/Tempnew - 1/Tempold) );
    w = w / sum(w);
    
    fprintf(1, 'Resampling acceptance rate ******************* %4.2f\n', resamp_accep/length(w));
    fprintf(1, 'Resampling rate of reduction ***************** %4.2f\n', ratenew);
    
end