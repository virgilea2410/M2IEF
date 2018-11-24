%% Computing the Jump Tests Statistics : 
%   - Test Statistics
%   - Pvalues 
%   - dataset = dataset of log-prices

function [all_scenar_haar_tests_stats, all_scenar_haar_pvalues, ...
            all_scenar_d4_tests_stats, all_scenar_d4_pvalues, ...
            all_scenar_s8_tests_stats, all_scenar_s8_pvalues, ...
            all_scenar_BNS_tests_stats, all_scenar_BNS_pvalues, ...
            all_scenar_JO_tests_stats, all_scenar_JO_pvalues, ...
            jump_sizes] ...
            = Jump_Test(dataset, compute_other_tests)

% Parameters
max_scenar = size(dataset, 2);
n_simul = size(dataset, 1);
%nu = round(1/sqrt(n_simul));
nu = 1;

% Test 12/06
%dataset = dataset - mean(dataset);
%dataset = exp(dataset);

% Price returns for BNS test statistic
price_returns_tmp = price2ret(exp(dataset));
price_returns = [zeros(1, size(dataset, 2));price_returns_tmp];

% Output Results
% Haar Results
all_haar_test_stats = zeros(n_simul, 1);
all_scenar_haar_tests_stats = zeros(n_simul, max_scenar);
all_haar_pvalues = zeros(n_simul, 1);
all_scenar_haar_pvalues = zeros(n_simul, max_scenar);
all_haar_jump_size = zeros(n_simul, 1);
all_scenar_haar_jump_size = zeros(n_simul, max_scenar);

% D4 Results
all_d4_test_stats = zeros(n_simul, 1);
all_scenar_d4_tests_stats = zeros(n_simul, max_scenar);
all_d4_pvalues = ones(n_simul, 1);
all_scenar_d4_pvalues = ones(n_simul, max_scenar);
all_d4_jump_size = zeros(n_simul, 1);
all_scenar_d4_jump_size = zeros(n_simul, max_scenar);

% S8 Results
all_s8_test_stats = zeros(n_simul, 1);
all_scenar_s8_tests_stats = zeros(n_simul, max_scenar);
all_s8_pvalues = ones(n_simul, 1);
all_scenar_s8_pvalues = ones(n_simul, max_scenar);
all_s8_jump_size = zeros(n_simul, 1);
all_scenar_s8_jump_size = zeros(n_simul, max_scenar);


% BSM (see Barndorff-Nielsen & Shepard (2006)) Results
all_BNS_test_stats = zeros(n_simul, 1);
all_scenar_BNS_tests_stats = zeros(n_simul, max_scenar);
all_BNS_pvalues = ones(n_simul, 1);
all_scenar_BNS_pvalues = ones(n_simul, max_scenar);
all_BNS_jump_size = zeros(n_simul, 1);
all_scenar_BNS_jump_size = zeros(n_simul, max_scenar);

% JO (see Jiang & Oomen (2008)) Results
all_JO_test_stats = zeros(n_simul, 1);
all_scenar_JO_tests_stats = zeros(n_simul, max_scenar);
all_JO_pvalues = ones(n_simul, 1);
all_scenar_JO_pvalues = ones(n_simul, max_scenar);
all_JO_jump_size = zeros(n_simul, 1);
all_scenar_JO_jump_size = zeros(n_simul, max_scenar);

for i_scenar=1:max_scenar
    for i_simul=3:n_simul
        % Test 14/06
        if i_simul == 400
            
        end
        % Wavelet Transformation of the vector for each wavelet filter
        modwt_haar = modwt(dataset(1:i_simul, i_scenar), 'haar');
        modwt_d4 = modwt(dataset(1:i_simul, i_scenar), 'db2');
        modwt_s8 = modwt(dataset(1:i_simul, i_scenar), 'sym4');

        if i_simul ~= n_simul
            Pplus = mean(exp(dataset(i_simul:i_simul+nu, i_scenar)));
            Pmoins = mean(exp(dataset(i_simul-nu:i_simul, i_scenar)));
        end
        
        % Test with Haar Filter Test Statistic
        P_haar = modwt_haar(1, end);
        integrated_vol_haar = 1/(i_simul-2) * sum(abs(modwt_haar(1,2:end-1)) .* ...
                                            abs(modwt_haar(1,1:end-2)));
        J_haar = (P_haar - mean(modwt_haar(1,1:end-1)))/sqrt(integrated_vol_haar);
        test_stat_haar = (J_haar)/sqrt(pi/2);
        pvalue_haar = 2*(1 - normcdf(abs(test_stat_haar)));
        all_haar_test_stats(i_simul, 1) = test_stat_haar;
        all_haar_pvalues(i_simul, 1) = pvalue_haar;
        all_haar_jump_size(i_simul, 1) = Pplus - Pmoins;
        
        % Test with D4 Filter Test Statistic
        P_d4 = modwt_d4(1, end);
        integrated_vol_d4 = 1/(i_simul-2) * sum(abs(modwt_d4(1,2:end-1)) .* ...
                                            abs(modwt_d4(1,1:end-2)));
        J_d4 = (P_d4 - mean(modwt_d4(1,1:end-1)))/sqrt(integrated_vol_d4);
        test_stat_d4 = (J_d4)/sqrt(3/2);
        pvalue_d4 = 2*(1 - normcdf(abs(test_stat_d4)));
        all_d4_test_stats(i_simul, 1) = test_stat_d4;
        all_d4_pvalues(i_simul, 1) = pvalue_d4;
        all_d4_jump_size(i_simul, 1) = Pplus - Pmoins;
        
        % Test with S8 Filter Test Statistic
        P_s8 = modwt_s8(1, end);
        integrated_vol_s8 = 1/(i_simul-2) * sum(abs(modwt_s8(1,2:end-1)) .* ...
                                            abs(modwt_s8(1,1:end-2)));
        J_s8 = (P_s8 - mean(modwt_s8(1,1:end-1)))/sqrt(integrated_vol_s8);
        test_stat_s8 = (J_s8)/1.18;
        pvalue_s8 = 2*(1 - normcdf(abs(test_stat_s8)));
        all_s8_test_stats(i_simul, 1) = test_stat_s8;
        all_s8_pvalues(i_simul, 1) = pvalue_s8;
        all_s8_jump_size(i_simul, 1) = Pplus - Pmoins;
        
        if compute_other_tests
            % Test with BNS Test Statistic
            mu1 = 2^(1/2) * gamma(1/2 * (1+1)) / gamma(1/2);
            test_stat_BNS_num = log(sum(price_returns(1:i_simul-1, i_scenar).^2)) - ...
                    log(1/mu1^2 * sum(abs(price_returns(1:i_simul-1, i_scenar)) .* ...
                    abs(price_returns(2:i_simul, i_scenar))));
            test_stat_BNS_denom = sqrt(0.6091 * max(1/i_simul, ...
                sum(abs(price_returns(1:i_simul-3, i_scenar)) .* abs(price_returns(2:i_simul-2, i_scenar)) ...
                .* abs(price_returns(3:i_simul-1, i_scenar)) .* abs(price_returns(4:i_simul, i_scenar))) ...
                / sum(abs((price_returns(1:i_simul-1, i_scenar))) .* abs(price_returns(2:i_simul, i_scenar)))^2));
            test_stat_BNS = test_stat_BNS_num / test_stat_BNS_denom;  
            pvalue_BNS = 2*(1 - normcdf(abs(test_stat_BNS)));    
            all_BNS_test_stats(i_simul, 1) = test_stat_BNS;
            all_BNS_pvalues(i_simul, 1) = pvalue_BNS;
            all_BNS_jump_size(i_simul, 1) = Pplus - Pmoins;

            % Test with JO Test Statistic
            SwVn = 2 * sum(price_returns(1:i_simul, i_scenar) - log(price_returns(1:i_simul, i_scenar) + 1));
            RVn = sum(log(price_returns(1:i_simul, i_scenar) + 1).^2);
            p = 4;
            mu6 = 2^(6/2) * gamma(1/2 * (6+1)) / gamma(1/2);
            mu6p = 2^((6/p)/2) * gamma(1/2 * ((6/p)+1)) / gamma(1/2);
            omega_SwV_4 = mu6/9 * n_simul^3 * mu6p^-p / (n_simul-p+1);
            tmp1 = 0;
            tmp2 = 1;
            for i=1:i_simul-p
                for k=1:p
                    tmp2 = tmp2 * abs(price_returns(i+k, i_scenar))^(6/p);
                end
                tmp1 = tmp1 + tmp2;
                tmp2 = 1;
            end
            omega_SwV_4 = omega_SwV_4 * tmp1;
            test_stat_JO = i_simul/sqrt(omega_SwV_4) * (SwVn - RVn);
            if abs(test_stat_JO) ~= Inf
                pvalue_JO = 2*(1 - normcdf(abs(test_stat_JO)));
            else
                pvalue_JO = -1;
                test_stat_JO = 0;
            end
            all_JO_test_stats(i_simul, 1) = test_stat_JO;
            all_JO_pvalues(i_simul, 1) = pvalue_JO;
            all_JO_jump_size(i_simul, 1) = Pplus - Pmoins;
        end
    end
    all_scenar_haar_tests_stats(:, i_scenar) = all_haar_test_stats;
    all_scenar_d4_tests_stats(:, i_scenar) = all_d4_test_stats;
    all_scenar_s8_tests_stats(:, i_scenar) = all_s8_test_stats;

    all_scenar_haar_pvalues(:, i_scenar) = all_haar_pvalues;
    all_scenar_d4_pvalues(:, i_scenar) = all_d4_pvalues;
    all_scenar_s8_pvalues(:, i_scenar) = all_s8_pvalues;
    
    all_scenar_haar_jump_size(:, i_scenar) = all_haar_jump_size;
    all_scenar_d4_jump_size(:, i_scenar) = all_d4_jump_size;
    all_scenar_s8_jump_size(:, i_scenar) = all_s8_jump_size;

    if compute_other_tests
        all_scenar_BNS_tests_stats(:, i_scenar) = all_BNS_test_stats;
        all_scenar_BNS_pvalues(:, i_scenar) = all_BNS_pvalues;
        all_scenar_BNS_jump_size(:, i_scenar) = all_BNS_jump_size;
        all_scenar_JO_tests_stats(:, i_scenar) = all_JO_test_stats;
        all_scenar_JO_pvalues(:, i_scenar) = all_JO_pvalues;
        all_scenar_JO_jump_size(:, i_scenar) = all_JO_jump_size;
    end
end

jump_sizes = all_haar_jump_size;

end