function [params, state_init]=params_base_condition()
params = struct();
params.R70 = 0.7; % fraction of total root in top 30 cm of soil
params.a = 600 / 365; % maximum nitrogen productivity gC/gN.day
params.b = 0.001; % the factor that decreases P_N = (a - bC_P)

%[Revision] Remove the logistic term. I made this value very large that effectively removes its effect
params.KNL = 75e8; % max leaf N
params.Vmax_rootCsupply = 1.5;% Root C supply rate constant gC leaf/gC root.d

% params.Vmax_exu_production = 0.0044*4 ; % maxmimum exudation rate constant d-1
params.Vmax_exu_production = 0.0019*4 ; % maxmimum exudation rate constant d-1 
%The global root exudate carbon flux doi: https://doi.org/10.1101/2024.02.01.578470 

params.tau_root = 0.4 / 365; % Root turnover rate constant d-1
params.Vmax_ex = 5.3/10; % Max. exudate uptake rate d-1
params.KM_ex = 3.3; % Half-saturation constant for exudate uptake kinetics gC/m2

params.Vmax_SOM = 0.15; % Max. SOM uptake rate d-1

params.KM_SOM = 6443; % Half saturation constant for SOM uptake rate gC/m2
params.tau_mic = 5 / 365; % Microbial turnover rate constant d-1

params.Vmax_IN_to_root = 3.7/15; % Max. root Inorganic N uptake rate d-1

params.KM_IN_to_root = 20; % Half saturation constant for root Inorganic N uptake rate
params.Vmax_IN_to_mic = 15; % Max. inorganic N supply rate to microbial pool d-1
params.KM_IN_to_mic = 1; % Half-saturation constant for inorganic N supply rate to microbial pool gN/m2
params.leaching_rate = 0.000; % inorganic N leaching_rate d-1
params.e_L = 0.7; % Leaf C assimilation efficiency
params.cue_R = 0.6; % Root carbon use efficiency
params.cue_M = 0.3; % microbial carbon use efficiency
params.CNE = 10000; % C-to-N ratio of exudates
params.CNM_clim = 10; % C limited C-to-N ratio of microbial pool

params.CNleaf_opt = 21; % Optimal C-to-N ratio of leaf

params.CNroot_opt = 50; %  C-to-N ratio of roots % range [10-100]
params.r_m = 0.005; % root maintenance respiration rate constant [0.001-0.01] (e.g. from Fan et al.: 0.72-4.32 mmol CO2-C g-1 DW d-1).
params.min_frac_leaf_growth = 0.2;
params.min_frac_root_growth = 0.2;
params.min_exudation =0.1;
params.max_exudation = 0.5;
% simulation parameters
% Define the parameters
init_root_shoot = 1.5; % range [0.3-1.8]
leafC = 50; % range [50-500]
rootC = init_root_shoot * leafC; % gC/m2

% Initialize the state struct
state_init = struct();

state_init.CNSOM0=10; % range [5-20];

state_init.leafC = leafC; % gC/m2
state_init.rootC = rootC; % gC/m2
state_init.exuC = 0.5; % gC/m2  % range [0.3-2] = [1.4-4.2] ugC/g soil
state_init.SOMC = 7200; % gC/m2  % range [1800-18000] == [0.5-5] %C
state_init.micC = 2.8; % gC/m2  % 0.04% of total SOC

% Calculate nitrogen states
state_init.leafN = state_init.leafC / params.CNleaf_opt; % gN/m2
state_init.rootN = state_init.rootC / params.CNroot_opt; % gN/m2
state_init.SOMN = state_init.SOMC / state_init.CNSOM0; % gN/m2
state_init.micN = state_init.micC / params.CNM_clim; % gN/m2
state_init.inorgN = 20; % gN/m2  % range 0.36 gN/m2 no fertilization to 20 gN/m2 with conventional fertilization
