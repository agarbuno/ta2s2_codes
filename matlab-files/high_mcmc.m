
%% Metropolis good start

close all

figure(1); clf; hold on
Nchain = 5000;
warmup = 500;
Ndim = 4;

for iter = 0:Ndim
    d = 10^iter;
    c = 2.38/sqrt(d);
    X = zeros(Nchain, d);

    x0 = randn(1,d);
    X(1,:) = x0;

    for i = 2:Nchain

        xnew = x0 + c^2 * randn(1, d);
        if rand < min(1, mvnpdf(xnew)/ mvnpdf(x0) ) 
            x0 = xnew;
        end
        X(i,:) = x0;

    end

    X = X(warmup:end,:);
    
    figure(1);
    h = histogram(sqrt(sum(X.^2,2))/sqrt(d));     
    h.Normalization = 'pdf';
    h.Parent.YScale = 'log';
    h.FaceAlpha = (iter+1)/5;
    h.FaceColor = 'blue';
    
    if d > 1 
        figure(4);
        plot(X(:,1), X(:,2), '.')
    end
    
    pause
end
hold off

%% Metropolis bad start

figure(2); clf; hold on

for iter = 0:Ndim
    d = 10^iter;
    c = 2.38/sqrt(d);
    X = zeros(Nchain, d);

    x0 = zeros(1,d);
    X(1,:) = x0;

    for i = 2:Nchain

        xnew = x0 + c^2 * randn(1, d);
        if rand < min(1, mvnpdf(xnew)/ mvnpdf(x0) ) 
            x0 = xnew;
        end
        X(i,:) = x0;

    end

    X = X(500:end,:);
    
    figure(2);
    h = histogram(sqrt(sum(X.^2,2))/sqrt(d));     
    h.Normalization = 'pdf';
    h.Parent.YScale = 'log';
    h.FaceAlpha = (iter+1)/5;
    h.FaceColor = 'blue';
    
    if d > 1 
        figure(4);
        plot(X(:,1), X(:,2), '.')
    end
    
    pause
end
hold off

%% Slice sampling

figure(3); clf; hold on

for iter = 0:Ndim

    d = 10^iter;
    c = 2.38/sqrt(d);
    X = zeros(Nchain, d);

    x0 = randn(1,d);
    X(1,:) = x0; pdf = mvnpdf(x0);

    for i = 1:(Nchain-1)

        % Defining the slice
        y = pdf * rand;

        % drawing first crumb
        ncrumbs = 1;
        crumb = x0 + c^2 * randn(1,d);
        xtry = crumb + c^2 * randn(1,d);
        pdftry = mvnpdf(xtry);

        while pdftry < y
            ncrumbs = ncrumbs + 1;
            crumbnew = x0 + c^2 * randn(1,d);
            crumb = 1/ncrumbs * ( (ncrumbs - 1) * crumb + crumbnew);

            xtry = crumb + (c/log(ncrumbs))^2 * randn(1,d);
            pdftry = mvnpdf(xtry);
        end

        X(i+1,:) = xtry;
        x0 = xtry; pdf = pdftry;

    end
    
    figure(3)
    h = histogram(sqrt(sum(X.^2,2))/sqrt(d));     
    h.Normalization = 'pdf';
    h.Parent.YScale = 'log';
    h.FaceAlpha = (iter+1)/5;
    h.FaceColor = 'blue';
    
    if d > 1 
        figure(4);
        plot(X(:,1), X(:,2), '.')
    end
    
    pause
end
hold off
    
%% Exact sampling

figure(5); clf; hold on
for i = 0:Ndim
    d = 10^i; X = randn( Nchain, d); 
    
    figure(5);
    h = histogram(sqrt(sum(X.^2,2))/sqrt(d));     
    h.Normalization = 'pdf';
    h.Parent.YScale = 'log';
    h.FaceAlpha = (1+i)/5;
    h.FaceColor = 'blue';
    
    if d > 1 
        figure(4);
        plot(X(:,1), X(:,2), '.')
    end
    
    pause
      
end
hold off
