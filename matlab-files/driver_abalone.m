% timevec = zeros(10,2);

for i = 9:10
    meanopt = 1;
    n = i * 100;
    execute_abalone
    timevec(i,1) = time; 
    abalone_residuals
    rmse(i,1) = sqrt(mean((mean_y - Ytry).^2));
    rmse_mix(i,1) = sqrt(mean((y_mix - Ytry).^2));
    save(sprintf('run_%03d.mat', i))
    clearvars -except timevec i rmse rmse_mix
end


% for i = 2:10
%     meanopt = 2;
%     n = i * 100;
%     execute_abalone
%     timevec(i,2) = time; 
%     abalone_residuals
%     rmse(i,2) = sqrt(mean((mean_y - Ytry).^2));
%     rmse_mix(i,2) = sqrt(mean((y_mix - Ytry).^2));
%     save(sprintf('run_linear_%03d.mat', i))
%     clearvars -except timevec i rmse rmse_mix
% end
