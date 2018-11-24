%% IV) Monte-Carlo Simulations
clc;
clear;
%% Under the null of no jumps

% Scenarios Parameters
mu_ranges = [0, 1];
vol_ranges = {'Constant', 'Stochastic'};

sample_size = 10000;
n_scenars = size(mu_ranges, 2) * size(vol_ranges, 2);
freq = 2;
%freq = 1;

%Alphas Errors
% Haar filter
all_haar_alphas_errors_5 = zeros(n_scenars, 1);
all_haar_alphas_errors_1 = zeros(n_scenars, 1);
% D4 Filter
all_d4_alphas_errors_5 = zeros(n_scenars, 1);
all_d4_alphas_errors_1 = zeros(n_scenars, 1);
% S8 Filter
all_s8_alphas_errors_5 = zeros(n_scenars, 1);
all_s8_alphas_errors_1 = zeros(n_scenars, 1);

TestSize = {'Variable_Name', 'Haar_Filter_1p100', 'D4_Filter_1p100', 'S8_Filter_1p100', 'Haar_Filter_5p100', 'D4_Filter_5p100', 'S8_Filter_5p100'; ...
            'Mean0_VolConst', [], [], [], [], [], []; ...
            'Mean0_VolSto', [], [], [], [], [], []; ...
            'Mean1_VolConst', [], [], [], [], [], []; ...
            'Mean1_VolSto', [], [], [], [], [], []};

% Simulates Stochastic Prices Without jumps
% & Computing Size
i_scenar = 1;
i_row = 2;
i_col = 2;
for i_mu = mu_ranges
    for i_vol = vol_ranges
        [stoch_price, stoch_log_price, stoch_vol, stoch_log_vol] ...
        = Simulate_Prices(sample_size, i_mu, i_vol, freq);

        % Computing the Jump Tests Statistics for simulated stochastic prices : 
        [all_scenar_haar_tests_stats, all_scenar_haar_pvalues, ...
            all_scenar_d4_tests_stats, all_scenar_d4_pvalues, ...
            all_scenar_s8_tests_stats, all_scenar_s8_pvalues] ...
            = Jump_Test(stoch_log_price, false);
        
        all_haar_alphas_errors_5(i_scenar, 1) = sum(all_scenar_haar_pvalues(:, 1) <= 0.05);
        all_d4_alphas_errors_5(i_scenar, 1) = sum(all_scenar_d4_pvalues(:, 1) <= 0.05);
        all_s8_alphas_errors_5(i_scenar, 1) = sum(all_scenar_s8_pvalues(:, 1) <= 0.05);
        
        all_haar_alphas_errors_1(i_scenar, 1) = sum(all_scenar_haar_pvalues(:, 1) <= 0.01);
        all_d4_alphas_errors_1(i_scenar, 1) = sum(all_scenar_d4_pvalues(:, 1) <= 0.01);
        all_s8_alphas_errors_1(i_scenar, 1) = sum(all_scenar_s8_pvalues(:, 1) <= 0.01);
        
        % Feeding the output table
        TestSize(i_row, i_col) = {all_haar_alphas_errors_1(i_scenar, 1)/sample_size};
        TestSize(i_row, i_col + 1) = {all_d4_alphas_errors_1(i_scenar, 1)/sample_size};
        TestSize(i_row, i_col + 2) = {all_s8_alphas_errors_1(i_scenar, 1)/sample_size};
        TestSize(i_row, i_col + 3) = {all_haar_alphas_errors_5(i_scenar, 1)/sample_size};
        TestSize(i_row, i_col + 4) = {all_d4_alphas_errors_5(i_scenar, 1)/sample_size};
        TestSize(i_row, i_col + 5) = {all_s8_alphas_errors_5(i_scenar, 1)/sample_size};
        
        % Plotting Test Stat For Haar Wavelet Filter
        figure()
        hold on;
        %[f, x] = ksdensity(all_scenar_haar_tests_stats(:, 1) * sqrt(pi/2));
        [f, x] = ksdensity(all_scenar_haar_tests_stats(:, 1));
        plot(x, f, 'r--o');
        plot(-4:0.1:4, normpdf(-4:0.1:4, 0, 1), 'g-');
        title('Test Statistic Density With Haar Wavelet Filter for a ' + string(i_vol) + ' Volatility, and Mean = ' +  string(i_mu));
        xlabel('Statistic Value');
        ylabel('Probability');
        legend({'Test Statistic', 'Standard Normal Distribution'});
        xlim([-5 5]);
        hold off;
        
        % Plotting Test Stat For D4 Wavelet Filter
        figure()
        hold on;
        %[f, x] = ksdensity(all_scenar_d4_tests_stats(:, 1) * sqrt(3/2));
        [f, x] = ksdensity(all_scenar_d4_tests_stats(:, 1));
        plot(x, f, 'r--o');
        plot(-4:0.1:4, normpdf(-4:0.1:4, 0, 1), 'g-');
        title('Test Statistic Density With D4 Wavelet Filter for a ' + string(i_vol) + ' Volatility, and Mean = ' +  string(i_mu));
        xlabel('Statistic Value');
        ylabel('Probability');
        legend({'Test Statistic', 'Standard Normal Distribution'});
        xlim([-5 5]);
        hold off;
        
        % Plotting Test Stat For S8 Wavelet Filter
        figure()
        hold on;
        %[f, x] = ksdensity(all_scenar_s8_tests_stats(:, 1) * sqrt(3/2));
        [f, x] = ksdensity(all_scenar_s8_tests_stats(:, 1));
        plot(x, f, 'r--o');
        plot(-4:0.1:4, normpdf(-4:0.1:4, 0, 1), 'g-');
        title('Test Statistic Density With S8 Wavelet Filter for a ' + string(i_vol) + ' Volatility, and Mean = ' +  string(i_mu));
        xlabel('Statistic Value');
        ylabel('Probability');
        legend({'Test Statistic', 'Standard Normal Distribution'});
        xlim([-5 5]);
        hold off;
        
        i_scenar = i_scenar + 1;
        i_col = 2;
        i_row = i_row + 1;
    end
end

TestSize = cell2table( ...
        TestSize(2:end, 2:end), ...
        'VariableNames', TestSize(1, 2:end), ...
        'RowNames', TestSize(2:end, 1));

% Displaying Size for each scnerios (Alpha, or Type I error)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Displaying Size of the test (alpha error)
disp('');
disp('------------------------- Size of the Test (Alpha Errors) ------------------------------');
disp(TestSize);
disp('');

%% Size and Power properties : under the alternative hypothesise of jump presence

% Computing Power (1 - Beta, or Type II error)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulate Stochastic Price with a jump
sample_size = 500;
n_simul = 50;
day_of_jump = round(0.8 * sample_size, 0);
range_freq = [3, 4];
range_size_jump = [0.1, 0.25, 0.5, 1, 3];

n_scenars = size(range_freq, 2) * size(range_size_jump, 2);

% Success Rate
% Haar Filter
all_haar_success_rate_5 = zeros(n_scenars, 1);
all_haar_success_rate_1 = zeros(n_scenars, 1);
% D4 Filter
all_d4_success_rate_5 = zeros(n_scenars, 1);
all_d4_success_rate_1 = zeros(n_scenars, 1);
% S8 Filter
all_s8_success_rate_5 = zeros(n_scenars, 1);
all_s8_success_rate_1 = zeros(n_scenars, 1);
% BNS Jump Test
all_BNS_success_rate_5 = zeros(n_scenars, 1);
all_BNS_success_rate_1 = zeros(n_scenars, 1);
% JO Jump Test
all_JO_success_rate_5 = zeros(n_scenars, 1);
all_JO_success_rate_1 = zeros(n_scenars, 1);

% Test Power
% Haar Filter
all_haar_test_power_5 = zeros(n_scenars, 1);
all_haar_test_power_1 = zeros(n_scenars, 1);
% D4 Filter
all_d4_test_power_5 = zeros(n_scenars, 1);
all_d4_test_power_1 = zeros(n_scenars, 1);
% S8 Filter
all_s8_test_power_5 = zeros(n_scenars, 1);
all_s8_test_power_1 = zeros(n_scenars, 1);
% BNS Filter
all_BNS_test_power_5 = zeros(n_scenars, 1);
all_BNS_test_power_1 = zeros(n_scenars, 1);
% JO Filter
all_JO_test_power_5 = zeros(n_scenars, 1);
all_JO_test_power_1 = zeros(n_scenars, 1);

% Boolean to compute power of the test
jump_is_detect_haar_5 = zeros(n_simul, 1);
jump_is_detect_d4_5 = zeros(n_simul, 1);
jump_is_detect_s8_5 = zeros(n_simul, 1);
jump_is_detect_BNS_5 = zeros(n_simul, 1);
jump_is_detect_JO_5 = zeros(n_simul, 1);

jump_is_detect_haar_1 = zeros(n_simul, 1);
jump_is_detect_d4_1 = zeros(n_simul, 1);
jump_is_detect_s8_1 = zeros(n_simul, 1);
jump_is_detect_BNS_1 = zeros(n_simul, 1);
jump_is_detect_JO_1 = zeros(n_simul, 1);

% Boolean to compute success rate of the test
no_spurious_jump_haar_5 = zeros(n_simul, 1);
no_spurious_jump_d4_5 = zeros(n_simul, 1);
no_spurious_jump_s8_5 = zeros(n_simul, 1);
no_spurious_jump_BNS_5 = zeros(n_simul, 1);
no_spurious_jump_JO_5 = zeros(n_simul, 1);

no_spurious_jump_haar_1 = zeros(n_simul, 1);
no_spurious_jump_d4_1 = zeros(n_simul, 1);
no_spurious_jump_s8_1 = zeros(n_simul, 1);
no_spurious_jump_BNS_1 = zeros(n_simul, 1);
no_spurious_jump_JO_1 = zeros(n_simul, 1);

TestPower = {'Variable_Name', 'S8_Filter', 'D4_Filter', 'Haar_Filter', 'BNS_Test', 'JO_Test'; ...
            'Sigma010_1H_1p100', [], [], [], [], []; ...
            'Sigma025_1H_1p100', [], [], [], [], []; ...
            'Sigma050_1H_1p100', [], [], [], [], []; ...
            'Sigma100_1H_1p100', [], [], [], [], []; ...
            'Sigma300_1H_1p100', [], [], [], [], []; ...
            'Sigma010_2H_1p100', [], [], [], [], []; ...
            'Sigma025_2H_1p100', [], [], [], [], []; ...
            'Sigma050_2H_1p100', [], [], [], [], []; ...
            'Sigma100_2H_1p100', [], [], [], [], []; ...
            'Sigma300_2H_1p100', [], [], [], [], []; ...
            'Sigma010_1H_5p100', [], [], [], [], []; ...
            'Sigma025_1H_5p100', [], [], [], [], []; ...
            'Sigma050_1H_5p100', [], [], [], [], []; ...
            'Sigma100_1H_5p100', [], [], [], [], []; ...
            'Sigma300_1H_5p100', [], [], [], [], []; ...
            'Sigma010_2H_5p100', [], [], [], [], []; ...
            'Sigma025_2H_5p100', [], [], [], [], []; ...
            'Sigma050_2H_5p100', [], [], [], [], []; ...
            'Sigma100_2H_5p100', [], [], [], [], []; ...
            'Sigma300_2H_5p100', [], [], [], [], []};
       
TestSuccessRate = {'Variable_Name', 'S8_Filter', 'D4_Filter', 'Haar_Filter', 'BNS_Test', 'JO_Test'; ...
            'Sigma010_1H_1p100', [], [], [], [], []; ...
            'Sigma025_1H_1p100', [], [], [], [], []; ...
            'Sigma050_1H_1p100', [], [], [], [], []; ...
            'Sigma100_1H_1p100', [], [], [], [], []; ...
            'Sigma300_1H_1p100', [], [], [], [], []; ...
            'Sigma010_2H_1p100', [], [], [], [], []; ...
            'Sigma025_2H_1p100', [], [], [], [], []; ...
            'Sigma050_2H_1p100', [], [], [], [], []; ...
            'Sigma100_2H_1p100', [], [], [], [], []; ...
            'Sigma300_2H_1p100', [], [], [], [], []; ...
            'Sigma010_1H_5p100', [], [], [], [], []; ...
            'Sigma025_1H_5p100', [], [], [], [], []; ...
            'Sigma050_1H_5p100', [], [], [], [], []; ...
            'Sigma100_1H_5p100', [], [], [], [], []; ...
            'Sigma300_1H_5p100', [], [], [], [], []; ...
            'Sigma010_2H_5p100', [], [], [], [], []; ...
            'Sigma025_2H_5p100', [], [], [], [], []; ...
            'Sigma050_2H_5p100', [], [], [], [], []; ...
            'Sigma100_2H_5p100', [], [], [], [], []; ...
            'Sigma300_2H_5p100', [], [], [], [], []};

i_scenar = 1;
i_row = 2;
i_col = 2;
for i_freq = range_freq
    for i_jump_size = range_size_jump
        for i_simul=1:n_simul
            % Simulates Stochastic Prices with a pre-defined jump
            [all_stoch_prices, all_stoch_log_prices, all_stoch_vols, all_stoch_log_vols] ...
                = Simulate_Prices_W_Jump(sample_size, 1, 'Constant', day_of_jump, i_jump_size, i_freq);

            % Computing the Jump Tests Statistics for simulated stochastic
            % prices with pre-defined jump : 
            [all_scenar_haar_tests_stats, all_scenar_haar_pvalues, ...
                all_scenar_d4_tests_stats, all_scenar_d4_pvalues, ...
                all_scenar_s8_tests_stats, all_scenar_s8_pvalues, ...
                all_scenar_BNS_tests_stats, all_scenar_BNS_pvalues, ...
                all_scenar_JO_tests_stats, all_scenar_JO_pvalues] ...
                = Jump_Test(all_stoch_log_prices, true);
    
            jump_is_detect_haar_5(i_simul, 1) = all_scenar_haar_pvalues(day_of_jump, 1) <= 0.05;
%             jump_is_detect_d4_5(i_simul, 1) = all_scenar_d4_pvalues(day_of_jump, 1) <= 0.05;
%             jump_is_detect_s8_5(i_simul, 1) = all_scenar_s8_pvalues(day_of_jump, 1) <= 0.05;
            jump_is_detect_d4_5(i_simul, 1) = all_scenar_d4_pvalues(day_of_jump+1, 1) <= 0.05;
            jump_is_detect_s8_5(i_simul, 1) = all_scenar_s8_pvalues(day_of_jump+1, 1) <= 0.05;
            jump_is_detect_BNS_5(i_simul, 1) = all_scenar_BNS_pvalues(day_of_jump, 1) <= 0.05;
            jump_is_detect_JO_5(i_simul, 1) = all_scenar_JO_pvalues(day_of_jump, 1) <= 0.05;
            
            jump_is_detect_haar_1(i_simul, 1) = all_scenar_haar_pvalues(day_of_jump, 1) <= 0.01;
%             jump_is_detect_d4_1(i_simul, 1) = all_scenar_d4_pvalues(day_of_jump, 1) <= 0.01;
%             jump_is_detect_s8_1(i_simul, 1) = all_scenar_s8_pvalues(day_of_jump, 1) <= 0.01;
            jump_is_detect_d4_1(i_simul, 1) = all_scenar_d4_pvalues(day_of_jump+1, 1) <= 0.01;
            jump_is_detect_s8_1(i_simul, 1) = all_scenar_s8_pvalues(day_of_jump+1, 1) <= 0.01;
            jump_is_detect_BNS_1(i_simul, 1) = all_scenar_BNS_pvalues(day_of_jump, 1) <= 0.01;
            jump_is_detect_JO_1(i_simul, 1) = all_scenar_JO_pvalues(day_of_jump, 1) <= 0.01;
            
            no_spurious_jump_haar_5(i_simul, 1) = (all_scenar_haar_pvalues(day_of_jump, 1) <= 0.05) * ...
                (sum(all_scenar_haar_pvalues(1:day_of_jump -1, 1) > 0.05) == 0) * ...
                (sum(all_scenar_haar_pvalues(day_of_jump+1:end, 1) > 0.05) == 0);
%             no_spurious_jump_d4_5(i_simul, 1) = (all_scenar_d4_pvalues(day_of_jump, 1) <= 0.05) * ...
%                 (sum(all_scenar_d4_pvalues(1:day_of_jump -1, 1) > 0.05) == 0) * ...
%                 (sum(all_scenar_d4_pvalues(day_of_jump+1:end, 1) > 0.05) == 0);
%             no_spurious_jump_s8_5(i_simul, 1) = (all_scenar_s8_pvalues(day_of_jump, 1) <= 0.05) * ...
%                 (sum(all_scenar_s8_pvalues(1:day_of_jump -1, 1) > 0.05) == 0) * ...
%                 (sum(all_scenar_s8_pvalues(day_of_jump+1:end, 1) > 0.05) == 0);
            no_spurious_jump_d4_5(i_simul, 1) = (all_scenar_d4_pvalues(day_of_jump+1, 1) <= 0.05) * ...
                (sum(all_scenar_d4_pvalues(1:day_of_jump, 1) > 0.05) == 0) * ...
                (sum(all_scenar_d4_pvalues(day_of_jump+2:end, 1) > 0.05) == 0);
            no_spurious_jump_s8_5(i_simul, 1) = (all_scenar_s8_pvalues(day_of_jump+1, 1) <= 0.05) * ...
                (sum(all_scenar_s8_pvalues(1:day_of_jump, 1) > 0.05) == 0) * ...
                (sum(all_scenar_s8_pvalues(day_of_jump+2:end, 1) > 0.05) == 0);
            no_spurious_jump_BNS_5(i_simul, 1) = (all_scenar_BNS_pvalues(day_of_jump, 1) <= 0.05) * ...
                (sum(all_scenar_BNS_pvalues(1:day_of_jump -1, 1) > 0.05) == 0) * ...
                (sum(all_scenar_BNS_pvalues(day_of_jump+1:end, 1) > 0.05) == 0);
            no_spurious_jump_JO_5(i_simul, 1) = (all_scenar_JO_pvalues(day_of_jump, 1) <= 0.05) * ...
                (sum(all_scenar_JO_pvalues(1:day_of_jump -1, 1) > 0.05) == 0) * ...
                (sum(all_scenar_JO_pvalues(day_of_jump+1:end, 1) > 0.05) == 0);

            no_spurious_jump_haar_1(i_simul, 1) = (all_scenar_haar_pvalues(day_of_jump, 1) <= 0.01) * ...
                (sum(all_scenar_haar_pvalues(1:day_of_jump -1, 1) > 0.01) == 0) * ...
                (sum(all_scenar_haar_pvalues(day_of_jump+1:end, 1) > 0.01) == 0);
%             no_spurious_jump_d4_1(i_simul, 1) = (all_scenar_d4_pvalues(day_of_jump, 1) <= 0.01) * ...
%                 (sum(all_scenar_d4_pvalues(1:day_of_jump -1, 1) > 0.01) == 0) * ...
%                 (sum(all_scenar_d4_pvalues(day_of_jump+1:end, 1) > 0.01) == 0);
%             no_spurious_jump_s8_1(i_simul, 1) = (all_scenar_s8_pvalues(day_of_jump, 1) <= 0.01) * ...
%                 (sum(all_scenar_s8_pvalues(1:day_of_jump -1, 1) > 0.01) == 0) * ...
%                 (sum(all_scenar_s8_pvalues(day_of_jump+1:end, 1) > 0.01) == 0);
            no_spurious_jump_d4_1(i_simul, 1) = (all_scenar_d4_pvalues(day_of_jump+1, 1) <= 0.01) * ...
                (sum(all_scenar_d4_pvalues(1:day_of_jump, 1) > 0.01) == 0) * ...
                (sum(all_scenar_d4_pvalues(day_of_jump+2:end, 1) > 0.01) == 0);
            no_spurious_jump_s8_1(i_simul, 1) = (all_scenar_s8_pvalues(day_of_jump+1, 1) <= 0.01) * ...
                (sum(all_scenar_s8_pvalues(1:day_of_jump, 1) > 0.01) == 0) * ...
                (sum(all_scenar_s8_pvalues(day_of_jump+2:end, 1) > 0.01) == 0);
            no_spurious_jump_BNS_1(i_simul, 1) = (all_scenar_BNS_pvalues(day_of_jump, 1) <= 0.01) * ...
                (sum(all_scenar_BNS_pvalues(1:day_of_jump -1, 1) > 0.01) == 0) * ...
                (sum(all_scenar_BNS_pvalues(day_of_jump+1:end, 1) > 0.01) == 0);
            no_spurious_jump_JO_1(i_simul, 1) = (all_scenar_JO_pvalues(day_of_jump, 1) <= 0.01) * ...
                (sum(all_scenar_JO_pvalues(1:day_of_jump -1, 1) > 0.01) == 0) * ...
                (sum(all_scenar_JO_pvalues(day_of_jump+1:end, 1) > 0.01) == 0);
        end
        all_haar_test_power_5(i_scenar, 1) = sum(jump_is_detect_haar_5 == 1)/n_simul;
        all_d4_test_power_5(i_scenar, 1) = sum(jump_is_detect_d4_5 == 1)/n_simul;
        all_s8_test_power_5(i_scenar, 1) = sum(jump_is_detect_s8_5 == 1)/n_simul;
        all_BNS_test_power_5(i_scenar, 1) = sum(jump_is_detect_BNS_5 == 1)/n_simul;
        all_JO_test_power_5(i_scenar, 1) = sum(jump_is_detect_JO_5 == 1)/n_simul;
        
        all_haar_test_power_1(i_scenar, 1) = sum(jump_is_detect_haar_1 == 1)/n_simul;
        all_d4_test_power_1(i_scenar, 1) = sum(jump_is_detect_d4_1 == 1)/n_simul;
        all_s8_test_power_1(i_scenar, 1) = sum(jump_is_detect_s8_1 == 1)/n_simul;
        all_BNS_test_power_1(i_scenar, 1) = sum(jump_is_detect_BNS_1 == 1)/n_simul;
        all_JO_test_power_1(i_scenar, 1) = sum(jump_is_detect_JO_1 == 1)/n_simul;
        
        all_haar_success_rate_5(i_scenar, 1) = sum(no_spurious_jump_haar_5 == 1)/n_simul;
        all_d4_success_rate_5(i_scenar, 1) = sum(no_spurious_jump_d4_5 == 1)/n_simul;
        all_s8_success_rate_5(i_scenar, 1) = sum(no_spurious_jump_s8_5 == 1)/n_simul;
        all_BNS_success_rate_5(i_scenar, 1) = sum(no_spurious_jump_BNS_5 == 1)/n_simul;
        all_JO_success_rate_5(i_scenar, 1) = sum(no_spurious_jump_JO_5 == 1)/n_simul;
        
        all_haar_success_rate_1(i_scenar, 1) = sum(no_spurious_jump_haar_1 == 1)/n_simul;
        all_d4_success_rate_1(i_scenar, 1) = sum(no_spurious_jump_d4_1 == 1)/n_simul;
        all_s8_success_rate_1(i_scenar, 1) = sum(no_spurious_jump_s8_1 == 1)/n_simul;
        all_BNS_success_rate_1(i_scenar, 1) = sum(no_spurious_jump_BNS_1 == 1)/n_simul;
        all_JO_success_rate_1(i_scenar, 1) = sum(no_spurious_jump_JO_1 == 1)/n_simul;
        
        TestPower(i_row, i_col) = {all_s8_test_power_1(i_scenar, 1)};
        TestPower(i_row, i_col + 1) = {all_d4_test_power_1(i_scenar, 1)};
        TestPower(i_row, i_col + 2) = {all_haar_test_power_1(i_scenar, 1)};
        TestPower(i_row, i_col + 3) = {all_BNS_test_power_1(i_scenar, 1)};
        TestPower(i_row, i_col + 4) = {all_JO_test_power_1(i_scenar, 1)};
        
        TestPower(i_row + 10, i_col) = {all_s8_test_power_5(i_scenar, 1)};
        TestPower(i_row + 10, i_col + 1) = {all_d4_test_power_5(i_scenar, 1)};
        TestPower(i_row + 10, i_col + 2) = {all_haar_test_power_5(i_scenar, 1)};
        TestPower(i_row + 10, i_col + 3) = {all_BNS_test_power_5(i_scenar, 1)};
        TestPower(i_row + 10, i_col + 4) = {all_JO_test_power_5(i_scenar, 1)};
        
        TestSuccessRate(i_row, i_col) = {all_s8_success_rate_1(i_scenar, 1)};
        TestSuccessRate(i_row, i_col + 1) = {all_d4_success_rate_1(i_scenar, 1)};
        TestSuccessRate(i_row, i_col + 2) = {all_haar_success_rate_1(i_scenar, 1)};
        TestSuccessRate(i_row, i_col + 3) = {all_BNS_success_rate_1(i_scenar, 1)};
        TestSuccessRate(i_row, i_col + 4) = {all_JO_success_rate_1(i_scenar, 1)};
        
        TestSuccessRate(i_row + 10, i_col) = {all_s8_success_rate_5(i_scenar, 1)};
        TestSuccessRate(i_row + 10, i_col + 1) = {all_d4_success_rate_5(i_scenar, 1)};
        TestSuccessRate(i_row + 10, i_col + 2) = {all_haar_success_rate_5(i_scenar, 1)};
        TestSuccessRate(i_row + 10, i_col + 3) = {all_BNS_success_rate_5(i_scenar, 1)};
        TestSuccessRate(i_row + 10, i_col + 4) = {all_JO_success_rate_5(i_scenar, 1)};
        
        i_scenar = i_scenar + 1;
        i_col = 2;
        i_row = i_row + 1;
    end
end

TestPower = cell2table( ...
        TestPower(2:end, 2:end), ...
        'VariableNames', TestPower(1, 2:end), ...
        'RowNames', TestPower(2:end, 1));

TestSuccessRate = cell2table( ...
        TestSuccessRate(2:end, 2:end), ...
        'VariableNames', TestSuccessRate(1, 2:end), ...
        'RowNames', TestSuccessRate(2:end, 1));    

% Displaying Power of the test (1-Beta error)
disp('---------------------------- Test Power -----------------------------');
disp(TestPower);
disp('---------------------------- Test Success Rate -----------------------------');
disp(TestSuccessRate);

%% V) Empirical Analysis for US Equity Markets

% Historical Data from 3 Stocks : General Electrics, IBM and Walmart
GE = xlsread('HistoPrices/GE_US.xlsx', 'Close');
IBM = xlsread('HistoPrices/IBM_US.xlsx', 'Close');
WMT = xlsread('HistoPrices/WMT_US.xlsx', 'Close');

% Dates and Time of historical datas
GE_1M_DATES = xlsread('HistoPrices/GE_US.xlsx', '1M', 'A1:A10000');
GE_5M_DATES = xlsread('HistoPrices/GE_US.xlsx', '5M', 'A1:A10000');
GE_15M_DATES = xlsread('HistoPrices/GE_US.xlsx', '15M', 'A1:A10000');

IBM_1M_DATES = xlsread('HistoPrices/IBM_US.xlsx', '1M', 'A1:A10000');
IBM_5M_DATES = xlsread('HistoPrices/IBM_US.xlsx', '5M', 'A1:A10000');
IBM_15M_DATES = xlsread('HistoPrices/IBM_US.xlsx', '15M', 'A1:A10000');

WMT_1M_DATES = xlsread('HistoPrices/WMT_US.xlsx', '1M', 'A1:A10000');
WMT_5M_DATES = xlsread('HistoPrices/WMT_US.xlsx', '5M', 'A1:A10000');
WMT_15M_DATES = xlsread('HistoPrices/WMT_US.xlsx', '15M', 'A1:A10000');

% Time of historical datas
GE_1M_TIME = datetime(datestr(GE_1M_DATES), 'Format', 'HH:mm:SS');
GE_5M_TIME = datetime(datestr(GE_5M_DATES), 'Format', 'HH:mm:SS');
GE_15M_TIME = datetime(datestr(GE_15M_DATES), 'Format', 'HH:mm:SS');

IBM_1M_TIME = datetime(datestr(IBM_1M_DATES), 'Format', 'HH:mm:SS');
IBM_5M_TIME = datetime(datestr(IBM_5M_DATES), 'Format', 'HH:mm:SS');
IBM_15M_TIME = datetime(datestr(IBM_15M_DATES), 'Format', 'HH:mm:SS');

WMT_1M_TIME = datetime(datestr(WMT_1M_DATES), 'Format', 'HH:mm:SS');
WMT_5M_TIME = datetime(datestr(WMT_5M_DATES), 'Format', 'HH:mm:SS');
WMT_15M_TIME = datetime(datestr(WMT_15M_DATES), 'Format', 'HH:mm:SS');

% Number of observation of datas
nobs_GE = size(GE, 1);
nobs_IBM = size(IBM, 1);
nobs_WMT = size(WMT, 1);

[GE_haar_tests_stats, GE_haar_pvalues, ...
GE_d4_tests_stats, GE_d4_pvalues, ...
GE_s8_tests_stats, GE_s8_pvalues, ...
~, ~, ~, ~, GE_jump_sizes] ...
= Jump_Test(log(GE), true);

[IBM_haar_tests_stats, IBM_haar_pvalues, ...
IBM_d4_tests_stats, IBM_d4_pvalues, ...
IBM_s8_tests_stats, IBM_s8_pvalues, ...
~, ~, ~, ~, IBM_jump_sizes] ...
= Jump_Test(log(IBM), true);

[WMT_haar_tests_stats, WMT_haar_pvalues, ...
WMT_d4_tests_stats, WMT_d4_pvalues, ...
WMT_s8_tests_stats, WMT_s8_pvalues, ...
~, ~, ~, ~, WMT_jump_sizes] ...
= Jump_Test(log(WMT), true);

% Total Jumps 
GE_haar_total_jump = sum(GE_haar_pvalues <= 0.05) ./ [nobs_GE/(60*24), nobs_GE/(12*24), nobs_GE/(4*24)];
IBM_haar_total_jump = sum(IBM_haar_pvalues <= 0.05) ./ [nobs_IBM/(60*24), nobs_IBM/(12*24), nobs_IBM/(4*24)];
WMT_haar_total_jump = sum(WMT_haar_pvalues <= 0.05) ./ [nobs_WMT/(60*24), nobs_WMT/(12*24), nobs_WMT/(4*24)];

GE_d4_total_jump = sum(GE_d4_pvalues <= 0.05) ./ [nobs_GE/(60*24), nobs_GE/(12*24), nobs_GE/(4*24)];
IBM_d4_total_jump = sum(IBM_d4_pvalues <= 0.05) ./ [nobs_IBM/(60*24), nobs_IBM/(12*24), nobs_IBM/(4*24)];
WMT_d4_total_jump = sum(WMT_d4_pvalues <= 0.05) ./ [nobs_WMT/(60*24), nobs_WMT/(12*24), nobs_WMT/(4*24)];

GE_s8_total_jump = sum(GE_s8_pvalues <= 0.05) ./ [nobs_GE/(60*24), nobs_GE/(12*24), nobs_GE/(4*24)];
IBM_s8_total_jump = sum(IBM_s8_pvalues <= 0.05) ./ [nobs_IBM/(60*24), nobs_IBM/(12*24), nobs_IBM/(4*24)];
WMT_s8_total_jump = sum(WMT_s8_pvalues <= 0.05) ./ [nobs_WMT/(60*24), nobs_WMT/(12*24), nobs_WMT/(4*24)];

% Total Positive Jumps
GE_haar_positive_jump = sum((GE_haar_pvalues <= 0.05) .* GE_jump_sizes > 0) ./ [nobs_GE/(60*24), nobs_GE/(12*24), nobs_GE/(4*24)];
IBM_haar_positive_jump = sum((IBM_haar_pvalues <= 0.05) .* IBM_jump_sizes > 0) ./ [nobs_IBM/(60*24), nobs_IBM/(12*24), nobs_IBM/(4*24)];
WMT_haar_positive_jump = sum((WMT_haar_pvalues <= 0.05) .* WMT_jump_sizes > 0) ./ [nobs_WMT/(60*24), nobs_WMT/(12*24), nobs_WMT/(4*24)];

GE_d4_positive_jump = sum((GE_d4_pvalues <= 0.05) .* GE_jump_sizes > 0) ./ [nobs_GE/(60*24), nobs_GE/(12*24), nobs_GE/(4*24)];
IBM_d4_positive_jump = sum((IBM_d4_pvalues <= 0.05) .* IBM_jump_sizes > 0) ./ [nobs_IBM/(60*24), nobs_IBM/(12*24), nobs_IBM/(4*24)];
WMT_d4_positive_jump = sum((WMT_d4_pvalues <= 0.05) .* WMT_jump_sizes > 0) ./ [nobs_WMT/(60*24), nobs_WMT/(12*24), nobs_WMT/(4*24)];

GE_s8_positive_jump = sum((GE_s8_pvalues <= 0.05) .* GE_jump_sizes > 0) ./ [nobs_GE/(60*24), nobs_GE/(12*24), nobs_GE/(4*24)];
IBM_s8_positive_jump = sum((IBM_s8_pvalues <= 0.05) .* IBM_jump_sizes > 0) ./ [nobs_IBM/(60*24), nobs_IBM/(12*24), nobs_IBM/(4*24)];
WMT_s8_positive_jump = sum((WMT_s8_pvalues <= 0.05) .* WMT_jump_sizes > 0) ./ [nobs_WMT/(60*24), nobs_WMT/(12*24), nobs_WMT/(4*24)];

% Total Negative Jumps
GE_haar_negative_jump = sum((GE_haar_pvalues <= 0.05) .* GE_jump_sizes < 0) ./ [nobs_GE/(60*24), nobs_GE/(12*24), nobs_GE/(4*24)];
IBM_haar_negative_jump = sum((IBM_haar_pvalues <= 0.05) .* IBM_jump_sizes < 0) ./ [nobs_IBM/(60*24), nobs_IBM/(12*24), nobs_IBM/(4*24)];
WMT_haar_negative_jump = sum((WMT_haar_pvalues <= 0.05) .* WMT_jump_sizes < 0) ./ [nobs_WMT/(60*24), nobs_WMT/(12*24), nobs_WMT/(4*24)];

GE_d4_negative_jump = sum((GE_d4_pvalues <= 0.05) .* GE_jump_sizes < 0) ./ [nobs_GE/(60*24), nobs_GE/(12*24), nobs_GE/(4*24)];
IBM_d4_negative_jump = sum((IBM_d4_pvalues <= 0.05) .* IBM_jump_sizes < 0) ./ [nobs_IBM/(60*24), nobs_IBM/(12*24), nobs_IBM/(4*24)];
WMT_d4_negative_jump = sum((WMT_d4_pvalues <= 0.05) .* WMT_jump_sizes < 0) ./ [nobs_WMT/(60*24), nobs_WMT/(12*24), nobs_WMT/(4*24)];

GE_s8_negative_jump = sum((GE_s8_pvalues <= 0.05) .* GE_jump_sizes < 0) ./ [nobs_GE/(60*24), nobs_GE/(12*24), nobs_GE/(4*24)];
IBM_s8_negative_jump = sum((IBM_s8_pvalues <= 0.05) .* IBM_jump_sizes < 0) ./ [nobs_IBM/(60*24), nobs_IBM/(12*24), nobs_IBM/(4*24)];
WMT_s8_negative_jump = sum((WMT_s8_pvalues <= 0.05) .* WMT_jump_sizes < 0) ./ [nobs_WMT/(60*24), nobs_WMT/(12*24), nobs_WMT/(4*24)];

% Replicate Jump Size matrix to match Pvalues matrix
GE_jump_sizes = repmat(GE_jump_sizes, 1, 3);
IBM_jump_sizes = repmat(IBM_jump_sizes, 1, 3);
WMT_jump_sizes = repmat(WMT_jump_sizes, 1, 3);

% Mean Amplitude of Jumps
GE_haar_mean_jump_size = mean(GE_jump_sizes .* (GE_haar_pvalues <= 0.05));
IBM_haar_mean_jump_size = mean(IBM_jump_sizes .* (IBM_haar_pvalues <= 0.05));
WMT_haar_mean_jump_size = mean(WMT_jump_sizes .* (WMT_haar_pvalues <= 0.05));

GE_d4_mean_jump_size = mean(GE_jump_sizes .* (GE_d4_pvalues <= 0.05));
IBM_d4_mean_jump_size = mean(IBM_jump_sizes .* (IBM_d4_pvalues <= 0.05));
WMT_d4_mean_jump_size = mean(WMT_jump_sizes .* (WMT_d4_pvalues <= 0.05));

GE_s8_mean_jump_size = mean(GE_jump_sizes .* (GE_s8_pvalues <= 0.05));
IBM_s8_mean_jump_size = mean(IBM_jump_sizes .* (IBM_s8_pvalues <= 0.05));
WMT_s8_mean_jump_size = mean(WMT_jump_sizes .* (WMT_s8_pvalues <= 0.05));

% Mean Amplitude of Positive Jumps
GE_haar_mean_pos_jump_size = mean(GE_jump_sizes .* (GE_haar_pvalues <= 0.05) .* (GE_jump_sizes > 0));
IBM_haar_mean_pos_jump_size = mean(IBM_jump_sizes .* (IBM_haar_pvalues <= 0.05) .* (IBM_jump_sizes > 0));
WMT_haar_mean_pos_jump_size = mean(WMT_jump_sizes .* (WMT_haar_pvalues <= 0.05) .* (WMT_jump_sizes > 0));

GE_d4_mean_pos_jump_size = mean(GE_jump_sizes .* (GE_d4_pvalues <= 0.05) .* (GE_jump_sizes > 0));
IBM_d4_mean_pos_jump_size = mean(IBM_jump_sizes .* (IBM_d4_pvalues <= 0.05) .* (IBM_jump_sizes > 0));
WMT_d4_mean_pos_jump_size = mean(WMT_jump_sizes .* (WMT_d4_pvalues <= 0.05) .* (WMT_jump_sizes > 0));

GE_s8_mean_pos_jump_size = mean(GE_jump_sizes .* (GE_s8_pvalues <= 0.05) .* (GE_jump_sizes > 0));
IBM_s8_mean_pos_jump_size = mean(IBM_jump_sizes .* (IBM_s8_pvalues <= 0.05) .* (IBM_jump_sizes > 0));
WMT_s8_mean_pos_jump_size = mean(WMT_jump_sizes .* (WMT_s8_pvalues <= 0.05) .* (WMT_jump_sizes > 0));

% Mean Amplitude of Negative Jumps
GE_haar_mean_neg_jump_size = mean(GE_jump_sizes .* (GE_haar_pvalues <= 0.05) .* (GE_jump_sizes < 0));
IBM_haar_mean_neg_jump_size = mean(IBM_jump_sizes .* (IBM_haar_pvalues <= 0.05) .* (IBM_jump_sizes < 0));
WMT_haar_mean_neg_jump_size = mean(WMT_jump_sizes .* (WMT_haar_pvalues <= 0.05) .* (WMT_jump_sizes < 0));

GE_d4_mean_neg_jump_size = mean(GE_jump_sizes .* (GE_d4_pvalues <= 0.05) .* (GE_jump_sizes < 0));
IBM_d4_mean_neg_jump_size = mean(IBM_jump_sizes .* (IBM_d4_pvalues <= 0.05) .* (IBM_jump_sizes < 0));
WMT_d4_mean_neg_jump_size = mean(WMT_jump_sizes .* (WMT_d4_pvalues <= 0.05) .* (WMT_jump_sizes < 0));

GE_s8_mean_neg_jump_size = mean(GE_jump_sizes .* (GE_s8_pvalues <= 0.05) .* (GE_jump_sizes < 0));
IBM_s8_mean_neg_jump_size = mean(IBM_jump_sizes .* (IBM_s8_pvalues <= 0.05) .* (IBM_jump_sizes < 0));
WMT_s8_mean_neg_jump_size = mean(WMT_jump_sizes .* (WMT_s8_pvalues <= 0.05) .* (WMT_jump_sizes < 0));

JumpDynamic = {'SamplingFrequency', 'Jump_Per_Day',     'Positive_Jump_Per_Day',        'Negative_Jump_Per_Day',                 'Mean_Jump_Size',               'Mean_Positive_Jump_Size',          'Mean_Negative_Jump_Size'; ...
                'GE_1M',    GE_haar_total_jump(1,1),     GE_haar_positive_jump(1,1),    GE_haar_negative_jump(1,1),     GE_haar_mean_jump_size(1, 1),   GE_haar_mean_pos_jump_size(1, 1),   GE_haar_mean_neg_jump_size(1, 1); ...
                'GE_5M',    GE_haar_total_jump(1,2),    GE_haar_positive_jump(1,2),     GE_haar_negative_jump(1,2),     GE_haar_mean_jump_size(1, 2),   GE_haar_mean_pos_jump_size(1, 2),   GE_haar_mean_neg_jump_size(1, 2); ...
                'GE_15M',   GE_haar_total_jump(1,3),    GE_haar_positive_jump(1,3),     GE_haar_negative_jump(1,3),     GE_haar_mean_jump_size(1, 3),   GE_haar_mean_pos_jump_size(1, 3),   GE_haar_mean_neg_jump_size(1, 3); ...
                'IBM_1M',   IBM_haar_total_jump(1,1),   IBM_haar_positive_jump(1,1),    IBM_haar_negative_jump(1,1),    IBM_haar_mean_jump_size(1, 1),  IBM_haar_mean_pos_jump_size(1, 1),  IBM_haar_mean_neg_jump_size(1, 1); ...
                'IBM_5M',   IBM_haar_total_jump(1,2),   IBM_haar_positive_jump(1,2),    IBM_haar_negative_jump(1,2),    IBM_haar_mean_jump_size(1, 2),  IBM_haar_mean_pos_jump_size(1, 2),  IBM_haar_mean_neg_jump_size(1, 2); ...
                'IBM_15M',  IBM_haar_total_jump(1,3),   IBM_haar_positive_jump(1,3),    IBM_haar_negative_jump(1,3),    IBM_haar_mean_jump_size(1, 3),  IBM_haar_mean_pos_jump_size(1, 3),  IBM_haar_mean_neg_jump_size(1, 3); ...
                'WMT_1M',   WMT_haar_total_jump(1,1),   WMT_haar_positive_jump(1,1),    WMT_haar_negative_jump(1,1),    WMT_haar_mean_jump_size(1, 1),  WMT_haar_mean_pos_jump_size(1, 1),  WMT_haar_mean_neg_jump_size(1, 1); ...
                'WMT_5M',   WMT_haar_total_jump(1,2),   WMT_haar_positive_jump(1,2),    WMT_haar_negative_jump(1,2),    WMT_haar_mean_jump_size(1, 2),  WMT_haar_mean_pos_jump_size(1, 2),  WMT_haar_mean_neg_jump_size(1, 2); ...
                'WMT_15M',  WMT_haar_total_jump(1,3),   WMT_haar_positive_jump(1,3),    WMT_haar_negative_jump(1,3),    WMT_haar_mean_jump_size(1, 3),  WMT_haar_mean_pos_jump_size(1, 3),  WMT_haar_mean_neg_jump_size(1, 3)};

JumpDynamic = cell2table( ...
        JumpDynamic(2:end, 2:end), ...
        'VariableNames', JumpDynamic(1, 2:end), ...
        'RowNames', JumpDynamic(2:end, 1)); 

disp('---------- Jump Dynamics for 3 major american stocks -----------')
disp(JumpDynamic)