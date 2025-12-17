The workflow:
1. Get input data which includes two types:
   (1) global datasets and undersampling in datasets by different depths and ocean basins ('Gptips_prerun_BAIT2_extractdata')
   (2) undersampling model output by using same depths and coords from IDP2025 ('Gptips_prerun_IDP_extractdata')

2. Run symbolic regression
   (1) Input data from 'Gptips_prerun_BAIT2_extractdata' should run in 'myconfig_BAIT2_fe3sol'
   (2) Input data from 'Gptips_prerun_IDP_extractdata' should run in 'myconfig_BAIT2_Obs_fe3sol'

3. Calculate EMD-SHAP
   (1) Equations discovered from global datasets (ideal scenario) should run in 'GPtips_postrun_SHAP_all_normalized'
   (2) Equations discovered from undersampling datasets should run in 'GPtips_postrun_SHAP_undersampling_normalized'
   Eventually then run in 'GPtips_postrun_EMD_all_normalized' to get EMD of SHAP between Millero equation and discovered equations.
