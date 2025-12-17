%% This script is to calculate SHAP of discovered equation from global dataset.
clear
close all
%% 1.Define symbolic variables and equations
syms tempis pH zsal;

% Define all 6 equations (replace with the actual equations)
% equation discovered from global datasets
eq1 = 0.02794529.*pH.^2 - 0.08564111.*pH - 0.10934083.*tempis + 0.01280994.*tempis.*pH.*exp(tempis) + 0.21783429;
eq2 = 0.03177317.*pH.*exp(tempis) - 0.11646833.*pH - 0.11565241.*tempis + 0.00312552.*pH.^4 + 0.22691366;
eq3 = 0.01927598.*pH.*exp(tempis) - 0.09770022.*pH - 0.12494497.*tempis + 0.01016874.*tempis.^3 + 0.02806899.*pH.^2 + 0.21436964;
eq4 = 0.01291578.*tempis.^3 + 0.03741823.*tempis.*pH - 0.13073923.*tempis - 0.01036830.*pH.^3 - 0.05697176.*pH + 0.22378780;
eq5 = 0.01421333.*tempis.^3 + 0.03274484.*tempis.*pH - 0.13106114.*tempis - 0.00637317.*pH.^3 + 0.01175895.*pH.^2 - 0.06006142.*pH + 0.21851058;
eq6 = 0.02514421.*tempis.*pH - 0.06726588.*pH - 0.13036370.*tempis + 0.01557831.*tempis.^3 + 0.01592579.*pH.^2 - 0.00503568.*tempis.*pH.*(pH - 0.19390000).*(1.64041272.*tempis - 1.64041272.*pH) + 0.22095899;

equations = {eq1, eq2, eq3, eq4, eq5, eq6};
n_equations = length(equations);

%% 2.Convert to Table-Compatible Functions
eq_functions = cell(1, n_equations);
for i = 1:n_equations
    % First create standard function
    func = matlabFunction(equations{i}, 'Vars', [tempis, zsal, pH]);

    % Wrap it to accept table input
    eq_functions{i} = @(tbl) func(tbl.tempis, tbl.zsal, tbl.pH);
end

%% 3. Load and Prepare Data
load("Gptips_CFe_BAIT2_openocean.mat"); % calculate SHAP value in the global dataset
load('Gptips_prerun_BAIT2_fe3sol_openocean_normalized.mat'); % replace here by different means and sigmas from different scenarios

% Create Latin Hypercube Sample for better coverage
n_samples = 1000;
rng(42); % For reproducibility

% Generate LHS sample in normalized space
lhs_sample = lhsdesign(n_samples, 3); % 3 variables
ranges = [min(BAIT2(:,[1,2,7])); max(BAIT2(:,[1,2,7]))]; % Physical ranges

% Scale LHS sample to physical ranges
X_phys = lhs_sample .* (ranges(2,:) - ranges(1,:)) + ranges(1,:);

% Normalize the sample
X_norm = (X_phys - mu_train) ./ sigma_train;

% Create input table
input_data = array2table(X_norm, 'VariableNames', {'tempis', 'zsal', 'pH'});

%% 4. Compute SHAP Values
shap_results = cell(n_equations, 1);
shap_metrics = struct('mean_abs_shap', [], 'mean_shap', [], 'std_shap', [], 'mean_contrib',[]);
shap_contrib = cell(n_equations,1);
shap_values = cell(n_equations,1);
threshold = 1e-4; % Values below this will be set to zero

for eq_idx = 1:n_equations
    fprintf('Processing equation %d/%d...\n', eq_idx, n_equations);

    % Create explainer
    rng(42);
    explainer = shapley(eq_functions{eq_idx}, input_data);

    % Compute SHAP values in parallel
    shap_results{eq_idx} = fit(explainer, input_data);

    % Calculate metrics
    vals = shap_results{eq_idx}.Shapley.Value; %3 x1000
    vals(abs(vals) < threshold) = 0; % Set small values to zero
    shap_contrib{eq_idx} = 100 .* abs(vals')./sum(abs(vals'),2); %1000 X3
    shap_metrics(eq_idx).mean_contrib = mean(shap_contrib{eq_idx},1);
    shap_metrics(eq_idx).mean_abs_shap = mean(abs(vals'), 1);
    shap_metrics(eq_idx).mean_shap = mean(vals', 1);
    shap_metrics(eq_idx).std_shap = std(vals', 0, 1);
    fprintf('Mean absolute SHAP: %.4f, %.4f, %.4f\n', shap_metrics(eq_idx).mean_abs_shap);

    %save for EMD
    shap_values{eq_idx} = vals;
    mean_shap(eq_idx,:) = shap_metrics(eq_idx).mean_shap;
    var_shap(eq_idx,:) = shap_metrics(eq_idx).std_shap;
end

% Save for EMD
contributions = shap_contrib;
save('/Users/chengwangwang/Documents/MATLAB/SHAP_test/shap_matlab_values_all_normalized.mat','n_equations','shap_values','mean_shap','var_shap','contributions');
