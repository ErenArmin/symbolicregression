function gp=myconfig_BAIT2_Obs_fe3sol(gp)

%% Run control
gp.runcontrol.pop_size = 250; % population of 250 models
gp.runcontrol.runs = 2; % perform 2 runs that are merged at the end
% gp.runcontrol.timeout = 60; % each run terminates after 60 seconds
% gp.runcontrol.parallel.auto = true;  % enable Parallel Computing if installed
gp.runcontrol.showBestInputs = true;
gp.runcontrol.showValBestInputs = true;

gp.fitness.terminate = true;           %true to enable early run termination on attaining a certain fitness value.
% gp.fitness.terminate_value = 0.001;      %terminate run early if this fitness value or better is achieved
gp.fitness.terminate_value = 0.015;      %terminate run early if this fitness value or better is achieved

% selection
gp.selection.tournament.size = 20; % usually under 10% of the population size
gp.selection.tournament.p_pareto = 1; % encourages less comples models
% gp.selection.tournament.p_pareto = 0.7; % encourages less comples models
% gp.selection.elite_fraction = 0.3; % approx. 1/3 models copied to next generation

% genes
gp.genes.max_genes = 6; % the max number of genes per model to 10
 
%% load data
load ("Gptips_CFe_BAIT2_openocean_IDP_DFe_IDP2025.mat");

[r,c] = size(BAIT2_obs) ; % r - rows, c - columns
P1 = 0.60 ; P2 = 0.80 ;  % P1 for training, (1-P2) for testing, P2-P1 for validation 
idx = randperm(r) ; % random number of rows, same length as 'r'
m = round(P1*r) ; n = round(P2*r) ;
BAIT2_train = BAIT2_obs(idx(1:m),:) ; 
BAIT2_val = BAIT2_obs(idx(m+1:n),:) ; 
BAIT2_test = BAIT2_obs(idx(n+1:end),:) ;

% Normalize using TRAINING stats only
% Compute ? and ? from training data
mu_train = mean(BAIT2_train(:, [1:2,7]));  % For T, S, pH
sigma_train = std(BAIT2_train(:, [1:2,7]));

BAIT2_train_normalized = BAIT2_train;
BAIT2_val_normalized = BAIT2_val;
BAIT2_test_normalized = BAIT2_test;

% Apply to all sets
BAIT2_train_normalized(:, [1:2,7]) = (BAIT2_train(:, [1:2,7]) - mu_train) ./ sigma_train;
BAIT2_val_normalized(:, [1:2,7]) = (BAIT2_val(:, [1:2,7]) - mu_train) ./ sigma_train;
BAIT2_test_normalized(:, [1:2,7]) = (BAIT2_test(:, [1:2,7]) - mu_train) ./ sigma_train;

save("/Users/chengwangwang/Documents/MATLAB/CFemodel/Gptips_prerun_BAIT2_fe3sol_openocean_IDP_DFe_IDP2025_normalized.mat",'BAIT2_train',...
    'BAIT2_val','BAIT2_test','BAIT2_train_normalized','BAIT2_val_normalized','BAIT2_test_normalized','mu_train','sigma_train');

%load train, test, validation data
gp.userdata.xtrain = BAIT2_train_normalized(:,[1:2,7]); % feature (X) T,S,pH
gp.userdata.ytrain = BAIT2_train(:,11); % output/target (Y) fe3sol

gp.userdata.xtest = BAIT2_test_normalized(:,[1:2,7]);
gp.userdata.ytest = BAIT2_test(:,11);

gp.userdata.xval = BAIT2_val_normalized(:,[1:2,7]);
gp.userdata.yval = BAIT2_val(:,11);

%enables hold out validation set
gp.userdata.user_fcn = @regressmulti_fitfun_validate;

%% fucntion definition
% give known variables aliases (this can include basic HTML markup)
gp.nodes.inputs.names = {'Temp','Sal','pH'};


%name for dataset
gp.userdata.name = 'PISCES-fe3sol-IDP-DFe-IDP2025-normalized'; 

% %define building block function nodes
% gp.nodes.functions.name = {'plus','minus','times','sqrt',...
%     'square','power'}; % trees should be built using which nodes

gp.nodes.functions.name = {'plus','minus','times','rdivide','cube','square',...
    'sqrt','exp','power','add3','mult3','log'}; % trees should be built using which nodes

% gp.nodes.functions.name = {'plus','minus','times','rdivide','cube','sqrt','abs',...
%     'square','exp','power','add3','mult3','log','negexp'}; % trees should be built using which nodes

end