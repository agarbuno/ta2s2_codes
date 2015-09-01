function [thetanew, Hnew, accep, ratenew, Ncrumbs] = ...
    slice_level(func, dim, c0, a, theta, ...
            N, w, rate, H, Temp, level, loc)

if N == 0   
    
    thetanew = []; 
    Hnew = []; 
    accep = []; 
    ratenew = [];
    
else
    
    accep = -1; nTOL = 100;

    % Adapting the spread
    ratenew = rate;
%     c = c0 * ratenew.^(level);
    c = c0;

    % Instead of selecting the most likely candidate to start the chain, use
    % the current loc level, this is for parallel chain
    index = loc;

    % Pre-allocate memory
    thetanew = zeros(N,size(theta,2));
    
    params = [theta];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  Using the covariance of the estimates %%%%
    if level == Inf 
        C = 5 * speye(size(params,2));
    else
        C = weightedcov(params, w);
    end
    % C = weightedcov(params, (1./H)/(sum(1./H)));
    try G = chol(C)';
    catch
        fprintf(1, '... Actually just use the variances ...\n');
        G = diag(std(params));
    end
%     G = diag(std(params));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [~, Ntheta] = size(theta);
    Hnew = zeros(N,1);

    % Initial state of the Markov Chain
    newparams = [theta(index,:)];
    
    % Extraction of the parameters
    thetanew(1,:) = newparams(1:Ntheta);
      
    Hnew(1,1) = feval(func, thetanew(1,:), a);

    % Compute the first slice 
    u = Hnew(1,1) + Temp * exprnd(1,1);

    for i = 1:(N+1)

        % Initializing the crumbs - avoiding the use of hyper rectangles
        ncrumbs = 1; 
        % Setting the current state of the chain
        m = [thetanew(i,:)];

        % Using only the appropriate crumbs for the annealing contrasted
        % with the current state of the Markov chain and the slice
        if (sum(H < u) > 1)  && (rand(1,1) < .5)
            % This is used to filter the crumbs inside the slice
            auxparams = params(H < u,:); auxw = w(H < u,:);
            crumb = auxparams(randsample(length(auxw), 1, true, auxw),:);
        else
            crumb = m + c0 * ( G * randn(Ntheta, 1))';
        end

        xiparams = crumb + c * ( G * randn(Ntheta, 1))';
         xitheta = xiparams(1:Ntheta);
           
        Hxi = feval(func, xitheta, a);
        
        % Resample until I have a good candidate using slice sampling
        while ((Hxi == Inf || Hxi == -Inf || imag(Hxi) ~= 0 || isnan(Hxi) ) || Hxi > u) && ...
                (ncrumbs < nTOL)
            ncrumbs = ncrumbs + 1;

            if (sum(H < u) > 1) && (rand(1,1) < .5)
                auxparams = params(H < u,:); auxw = w(H < u,:);
                crumbnew = auxparams(randsample(length(auxw), 1, true, auxw),:);
            else
                crumbnew = m + c0 * ( G * randn(Ntheta, 1))';
            end

            crumb = 1/ncrumbs * ( (ncrumbs - 1) * crumb + crumbnew);
            alphacrumb = 1 - 1/ncrumbs; 

            crumb = alphacrumb * m + (1 - alphacrumb) * crumb;

            xiparams = crumb + (c/ncrumbs)  * ( G * randn(Ntheta, 1))';
             xitheta = xiparams(1:Ntheta);
               
            Hxi = feval(func, xitheta, a);
            
        end

        Ncrumbs(i+1,1) = ncrumbs;

        if ncrumbs >= nTOL 
            thetanew(i+1, :) = thetanew(i, :);
              Hnew(i+1,1) = Hnew(i, :);
        else 
            thetanew(i+1, :) = xitheta;
              Hnew(i+1,1) = Hxi;

                 accep = accep + 1;
        end

        u = Hnew(i+1,1) + Temp * exprnd(1,1);

    end
    
    Ncrumbs = Ncrumbs(2:(N+1),:);
    thetanew = thetanew(2:(N+1),:);
    Hnew = Hnew(2:(N+1),:);

end