%% This script is to calculate the EMD of SHAP between Millero equation and other discovered equations.
close all
clear
%% Compute EMD of SHAP between Millero equations and other discovered equations

% can replace here by different SHAP files from different undersampling scenarioes
load('/Users/chengwangwang/Documents/MATLAB/SHAP_test/shap_matlab_values_all_normalized.mat'); % including numbers of equations, shap, its var and mean and features' contribution
load('shap_matlab_values_eq0_normalized.mat'); % eq0 is Millero

% Convert each 3x1000 matrix to 1000x3 in the cell array
for eq_idx = 1:n_equations
    shap_values{eq_idx} = shap_values{eq_idx}';  % Transpose each matrix
end
shap_values_eq0{1} = shap_values_eq0{1}';

% EMD results (n_equations x 3 features)
emd_results = zeros(n_equations, 3);

% For each equation (eq1 to eq14)
for eq_idx = 1:n_equations
    % For each feature (tempis, pH, zsal)
    for feature_idx = 1:3
        % Get SHAP values for current feature
        shap_eq0 = shap_values_eq0{1}(:, feature_idx);
        shap_eq = shap_values{eq_idx}(:, feature_idx);

        % Remove NaN/Inf (critical for stability)
        valid_mask = ~isnan(shap_eq0) & ~isinf(shap_eq0) & ~isnan(shap_eq) & ~isinf(shap_eq);
        shap_eq0 = shap_eq0(valid_mask);
        shap_eq = shap_eq(valid_mask);

        % Compute Wasserstein Distance (EMD) directly on raw values
        emd_results(eq_idx, feature_idx) = GPtips_postrun_wasserstein_distance(shap_eq0, shap_eq,50);
    end
end

%% Cleaner Implementation of Weighted Sums Calculation
% Initialize matrices
weighted_sums = zeros(n_equations, 3);      % Individual variable impacts
combined_weighted_sums = zeros(n_equations, 1); % Combined impacts

for eq_idx = 1:n_equations
    mean_contrib = mean(contributions{eq_idx}, 1); % [tempis, pH, zsal] contributions
    emd = emd_results(eq_idx,:);

    % % Calculate contrib × SHAP in one vectorized operation
    weighted_sums(eq_idx,:) = emd;
    combined_weighted_sums(eq_idx) = sum(emd);

end

%% Display Results Cleanly
variables = {'tempis', 'zsal', 'pH'};

results_table = array2table([weighted_sums, combined_weighted_sums], ...
    'VariableNames', [variables, {'Combined_EMD'}], ...
    'RowNames', arrayfun(@(x) sprintf('Eq%d', x), 1:n_equations, 'UniformOutput', false));

disp('EMD of SHAP:');
disp(results_table);
