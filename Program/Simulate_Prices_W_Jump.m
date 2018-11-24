%% Simulating Stochastic Prices With a Jump at specified observation :
%   - With Specified Drift
%   - With Stochastic or Constant Volatiliy
%   - With specified Jump localisation
%   - With specified Jump size

function [stoch_price, stoch_log_price, ...
            stoch_var, stoch_log_var] = ...
            Simulate_Prices_W_Jump(n_simul, price_mean, vol_type, t_jump, size_jump, frequency)
        
% Choosen Parameters
dt = 0;
if frequency == 1
    % 1 minute tick
    dt = 1/(60*24*365);
elseif frequency == 2
    % 15 minutes frequency
    dt = 1/(4*24*365);
elseif frequency == 3
    % 1 hour frequency
    dt = 1/(24*365);
elseif frequency == 4
    % 2 hour frequency
    dt = 1/(12*365);
end

% Given Parameters
kappa = -0.1;
delta = 0.25;
theta = -6.802;
rho = 0.5;

hist_vol = sqrt(exp(theta));
hist_price = 100;

% Output Vectors & Matrixs
stoch_log_var = zeros(n_simul, 1);
stoch_log_price = zeros(n_simul, 1);
stoch_var = zeros(n_simul, 1);
stoch_price = zeros(n_simul, 1);

% Initialisazing P0 and V0 (Start value for Price & Volatility)
stoch_log_price(1,1) = log(hist_price);
stoch_log_var(1,1) = log(hist_vol^2);
stoch_price(1,1) = hist_price;
stoch_var(1,1) = hist_vol^2;

% Simulating Stochastic Prices
if strcmp(vol_type, 'Constant')
    for i_simul=2:n_simul
        if i_simul ~= t_jump
            dWS = randn;
%         stoch_log_price(i_simul, 1) = stoch_log_price(i_simul - 1, 1) + ...
%             (price_mean - 1/2 * hist_vol^2) * dt + hist_vol * sqrt(dt) * dWS;
            stoch_log_price(i_simul, 1) = stoch_log_price(i_simul - 1, 1) + ...
                (price_mean - 1/2 * hist_vol) * dt + sqrt(hist_vol) * sqrt(dt) * dWS;
            stoch_price(i_simul, 1) = exp(stoch_log_price(i_simul, 1));                
        else
            dWS = randn;
%         stoch_log_price(i_simul, 1) = stoch_log_price(i_simul - 1, 1) + ...
%             (price_mean - 1/2 * hist_vol^2) * dt + hist_vol * sqrt(dt) * dWS;
        stoch_log_price(i_simul, 1) = stoch_log_price(i_simul - 1, 1) + ...
            (price_mean - 1/2 * hist_vol) * dt + sqrt(hist_vol) * sqrt(dt) * dWS;
            stoch_price(i_simul, 1) = exp(stoch_log_price(i_simul, 1)) + size_jump * hist_vol;
            stoch_price(i_simul, 1) = exp(stoch_log_price(i_simul, 1)) + size_jump * std(exp(stoch_log_price(1:i_simul-1, 1)));
            stoch_log_price(i_simul, 1) = log(stoch_price(i_simul, 1));
        end
    end
elseif strcmp(vol_type, 'Stochastic')
    for i_simul=2:n_simul
        if i_simul ~= t_jump
            dWS = randn;
            dWV = randn;
                
            stoch_log_var(i_simul, 1) = stoch_log_var(i_simul - 1, 1) + kappa * ...
                (theta - stoch_log_var(i_simul - 1, 1)) * dt + delta * sqrt(stoch_var(i_simul - 1, 1)) * sqrt(dt) * dWV;
            stoch_var(i_simul, 1) = exp(stoch_log_var(i_simul, 1));
                
            stoch_log_price(i_simul, 1) = stoch_log_price(i_simul - 1, 1) + ...
                (price_mean - 1/2 * max(sqrt(stoch_var(i_simul, 1)), 0)) * dt + ...
                sqrt(max(sqrt(stoch_var(i_simul, 1)), 0)) * sqrt(dt) * dWS;
            stoch_price(i_simul, 1) = exp(stoch_log_price(i_simul, 1));
        else
            dWS = randn;
            dWV = randn;
                
            stoch_log_var(i_simul, 1) = stoch_log_var(i_simul - 1, 1) + kappa * ...
                (theta - stoch_log_var(i_simul - 1, 1)) * dt + delta * sqrt(stoch_var(i_simul - 1, 1)) * sqrt(dt) * dWV;
            stoch_var(i_simul, 1) = exp(stoch_log_var(i_simul, 1));
                
            stoch_log_price(i_simul, 1) = stoch_log_price(i_simul - 1, 1) + ...
                (price_mean - 1/2 * max(sqrt(stoch_var(i_simul, 1)), 0)) * dt + ...
                sqrt(max(sqrt(stoch_var(i_simul, 1), 0))) * sqrt(dt) * dWS;
            stoch_price(i_simul, 1) = exp(stoch_log_price(i_simul, 1)) + ...
                size_jump * max(stoch_var(i_simul, 1), 0);
            stoch_log_price(i_simul, 1) = log(stoch_price(i_simul, 1));
        end
    end
end
end
