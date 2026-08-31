function [df, sol, litter, ocp] = solve_ocp_nested(terminal_time, state_init, params, num_tstep, tdf)
% SOLVE_OCP_NESTED - Optimal control with nested (coupled) control structure
%
% This version implements Option 1: Hierarchical control where exudation is
% a fraction of net root C assimilation, ensuring biological coupling between
% root C supply (S) and exudation (E).
%
% Control structure:
%   u(1) = rootCSupply [gC/m²/d] - C allocated from leaves to roots
%   u(2) = exu_fraction [-] - Fraction of net C supplied to root allocated to exudation
%   
%   Then: root_exu = exu_fraction * Net_rootC_assimilation
%
% This prevents E and S from competing artificially and eliminates bang-bang
% control behavior.

% Create the Yop system
sys = YopSystem( ...
    'states', 9, ...
    'controls', 2, ...
    'model', @(t, x, u)plant_soil_microbe_nested(t, x, u, params) ...
    );

% Symbolic variables
time = sys.t;
litter = sys.y.litter;
x = sys.x;
u = sys.u;

% Formulate optimal control problem
ocp = YopOcp();
% ocp.max({timeIntegral(litter.leaf_growth_rate) });

% ocp.max({timeIntegral(litter.leaf_growth_rate) '+' timeIntegral((diff(u(2)))^2)* 1e-6 });
% First term: Integral of leaf growth rate
term1 = timeIntegral(litter.leaf_growth_rate);
ocp.max({ term1 })

ocp.st( ...
    'systems', sys, ...
    ... % State bounds
    {0, '<=', litter.leafC}, ...
    {0, '<=', litter.rootC}, ...
    {0, '<=', litter.exuC}, ...
    {0, '<=', litter.SOMC}, ...
    {0, '<=', litter.micC}, ...
    {0, '<=', litter.rootN}, ...
    {0, '<=', litter.SOMN}, ...
    {0, '<=', litter.micN}, ...
    {0.01, '<=', litter.inorgN}, ...
    ... % Control bounds
    {0, '<=', litter.CNM, '<=', 20}, ...
    {0, '<=', litter.leaf_growth_rate}, ...
    {0, '<=', litter.root_exu}, ...
    {0, '<=', litter.rootCSupply}, ...
    ... % NEW: Bound exu_fraction instead of absolute exudation
    {params.min_exudation, '<=', u(2), '<=', params.max_exudation}, ...  % exu_fraction ∈ [0.1, 0.5]
    {0, '<=', litter.rootN_uptake_norm, '<=', 1}, ...
    {params.min_frac_leaf_growth, '<=', litter.frac_leaf_growth}, ...
    {params.min_frac_root_growth, '<=', litter.frac_root_growth}, ...
    ... % Initial conditions
    {0, '==', t_0(time)}, ...
    {state_init.leafC, '==', t_0(litter.leafC)}, ...
    {state_init.rootC, '==', t_0(litter.rootC)}, ...
    {state_init.exuC, '==', t_0(litter.exuC)}, ...
    {state_init.SOMC, '==', t_0(litter.SOMC)}, ...
    {state_init.micC, '==', t_0(litter.micC)}, ...
    {state_init.rootN, '==', t_0(litter.rootN)}, ...
    {state_init.SOMN, '==', t_0(litter.SOMN)}, ...
    {state_init.micN, '==', t_0(litter.micN)}, ...
    {state_init.inorgN, '==', t_0(litter.inorgN)}, ...
    ... % Terminal conditions
    {terminal_time, '==', t_f(time)} ...
    );

% Solving the OCP
if ~isempty(tdf)
    % Warm start: interpolate previous solution
    % Need to convert root_exu back to exu_fraction for initial guess
    state = tdf{:, ["leafC", "rootC", "exuC", "SOMC", "micC", "rootN", "SOMN", "micN", "inorgN", "rootCSupply"]};
    
    % Estimate exu_fraction from previous run
    if ismember('root_exu', tdf.Properties.VariableNames) && ismember('Net_rootC_assimilation', tdf.Properties.VariableNames)
        exu_frac_guess = tdf.root_exu ./ max(tdf.Net_rootC_assimilation, 0.01);
        exu_frac_guess = max(params.min_exudation, min(0.5, exu_frac_guess));  % Clip to bounds
    else
        exu_frac_guess = ones(height(tdf), 1) * 0.25;  % Default guess
    end
    
    state = [state, array2table(exu_frac_guess, 'VariableNames', {'exu_fraction'})];
    state_t = tdf.time;
    new_t = linspace(0, terminal_time, num_tstep + 1);
    new_state = interp1(state_t, state{:, :}, new_t, "pchip")';
    sz = size(new_state);
    new_state = new_state + new_state .* rand(sz) * 0.1;
    InitialGuess = [new_t; new_state];
    
    w0 = YopInitialGuess(...
        'signals', [time; x; u], ...
        'signalValues', InitialGuess ...
        );
    
    tic
    sol = ocp.solve( ...
        'controlIntervals', num_tstep, ...
        'collocationPoints', 'legendre', ...
        'polynomialDegree', 5, ...
        'initialGuess', w0, ...
        'ipopt', struct('max_iter', 2000, "print_level", 3, "max_cpu_time", 500,...
        'nlp_scaling_method','gradient-based') ...
        );
    toc
else
    % Cold start
    tic
    sol = ocp.solve( ...
        'controlIntervals', num_tstep, ...
        'collocationPoints', 'legendre', ...
        'polynomialDegree', 3, ...
        'ipopt', struct('max_iter', 2000, "print_level", 3, "max_cpu_time", 500, ...
        'tol', 1e-6, 'warm_start_init_point', 'yes', 'output_file', 'out.txt') ...
        );
    toc
end

% Extract solution
df = array2table(sol.signal(time)', 'VariableNames', {'time'});
litter_cell = struct2cell(litter);
fields = fieldnames(litter);
for i = 1:length(fields)
    varName = fields{i};
    df.(varName) = sol.signal(litter_cell{i})';
end
df.CNE = ones(height(df), 1) * params.CNE;

% Add exu_fraction control to output
df.exu_fraction = sol.signal(u(2))';

    function [dx, y] = plant_soil_microbe_nested(t, state, u, params)
        % PLANT_SOIL_MICROBE_NESTED - ODE system with nested control structure
        %
        % Key change: u(2) is now exu_fraction (dimensionless) instead of
        % root_exu (absolute rate). Exudation is computed as:
        %   root_exu = exu_fraction * Net_rootC_assimilation
        
        % Define the state variables
        fCL = state(1);
        fCR = state(2);
        fCE = state(3);
        fSOMC = state(4);
        fCM = state(5);
        
        fNL = fCL / params.CNleaf_opt;
        fNR = state(6);
        fNSOM = state(7);
        fNM = state(8);
        fNinorg = state(9);
        
        % Extract controls
        rootCSupply = u(1);      % Root C supply [gC/m²/d]
        exu_fraction = u(2);     % Exudation fraction [-]
        
        % Calculate various ratios
        CNR = fCR / fNR;
        CNM = fCM / fNM;
        CNSOM = fSOMC / fNSOM;
        CNleaf = fCL / fNL;
        
        % Calculate C assimilation rate A
        A = fNL * (1 - fNL / params.KNL) * (params.a - params.b * fCL);
        
        % Calculate leaf respiration
        leaf_resp = (1 - params.e_L) * A;
        
        % Calculate root turnover
        root_turnover = params.tau_root * fCR;
        
        % Calculate exudate uptake
        exu_uptake = (params.Vmax_ex * fCE * fCM) / (params.KM_ex + fCM + fCE);
        
        % Calculate SOM uptake
        SOM_uptake = (params.Vmax_SOM * fSOMC * fCM) / (params.KM_SOM + fCM + fSOMC);
        
        % Calculate microbial turnover
        micro_turnover = params.tau_mic * fCM;
        
        % Calculate microbial respiration
        micro_resp = (1 - params.cue_M) * (exu_uptake + SOM_uptake);
        
        % Calculate potential root N uptake
        potential_root_N_uptake = (params.Vmax_IN_to_root * fNR * fNinorg) / (params.KM_IN_to_root + fNinorg);
        
        % Calculate root respiration
        root_resp = (1 - params.cue_R) * rootCSupply;
        root_maintenance_resp = params.r_m * fCR;
        
        % Calculate leaf N uptake
        leaf_N_uptake = (A - leaf_resp - rootCSupply) / params.CNleaf_opt;
        
        % Net root C available for allocation
        Net_rootC_assimilation = rootCSupply - root_resp;
        
        % KEY CHANGE: Compute exudation from fraction
        root_exu = exu_fraction * Net_rootC_assimilation;
        
        % Root growth is remainder after exudation and maintenance
        root_growth_rate = Net_rootC_assimilation - root_exu - root_maintenance_resp;
        
        % Calculate phiC (microbial C-N coupling)
        phiC = SOM_uptake / CNSOM + exu_uptake / params.CNE - ...
            params.cue_M * (exu_uptake + SOM_uptake) / CNM;
        
        % Calculate inorganic N input to microbial pool
        IN_inorg_to_mic = params.Vmax_IN_to_mic * fNinorg / (params.KM_IN_to_mic + fNinorg);
        
        % Calculate N leaching
        Nleaching = params.leaching_rate * fNinorg;
        
        % N mineralization (conditional on phiC sign)
        Nmin = if_else(phiC > 0, phiC, 0);
        
        % N immobilization (conditional logic)
        Imm = if_else(phiC > 0, 0, if_else(-phiC > IN_inorg_to_mic, IN_inorg_to_mic, -phiC));
        
        % Net N mineralization
        phiN = Nmin - Imm;
        
        % Leaf C budget
        Net_leafC_assimilation = A - leaf_resp;
        leaf_growth_rate = Net_leafC_assimilation - rootCSupply;
        frac_Anet_allocated_toleaf_growth = leaf_growth_rate / Net_leafC_assimilation;
        
        % Root allocation fractions
        root_CUE = Net_rootC_assimilation / rootCSupply;
        frac_S_allocated_to_root_growth = root_growth_rate / Net_rootC_assimilation;
        
        % Root N uptake
        root_N_uptake = root_growth_rate / params.CNroot_opt + leaf_N_uptake + root_exu / params.CNE;
        
        % State derivatives
        dleafCdt = leaf_growth_rate;
        drootCdt = root_growth_rate - root_turnover;
        dexuCdt = root_exu - exu_uptake;
        dSOMCdt = -SOM_uptake + micro_turnover + root_turnover;
        dmicCdt = exu_uptake + SOM_uptake - micro_resp - micro_turnover;
        
        drootNdt = root_N_uptake - root_exu / params.CNE - root_turnover / params.CNroot_opt - leaf_N_uptake;
        dSOMNdt = -SOM_uptake / CNSOM + micro_turnover / CNM + root_turnover / params.CNroot_opt;
        dmicNdt = exu_uptake / params.CNE + SOM_uptake / CNSOM - micro_turnover / CNM - phiN;
        dInorgNdt = 0 + phiN - root_N_uptake - Nleaching;
        
        dx = [dleafCdt; drootCdt; dexuCdt; dSOMCdt; dmicCdt; ...
            drootNdt; dSOMNdt; dmicNdt; dInorgNdt];
        
        % Output variables
        y.litter.leafC = fCL;
        y.litter.rootC = fCR;
        y.litter.exuC = fCE;
        y.litter.SOMC = fSOMC;
        y.litter.micC = fCM;
        
        y.litter.leafN = fNL;
        y.litter.rootN = fNR;
        y.litter.Nexu = fCE / params.CNE;
        y.litter.SOMN = fNSOM;
        y.litter.micN = fNM;
        y.litter.inorgN = fNinorg;

        y.litter.S_norm = rootCSupply / Net_leafC_assimilation;
        y.litter.frac_leaf_growth = frac_Anet_allocated_toleaf_growth;
        y.litter.root_CUE = root_CUE;
        y.litter.frac_root_growth = frac_S_allocated_to_root_growth;
        
        % Ratios
        y.litter.CNR = CNR;
        y.litter.CNM = CNM;
        y.litter.CNSOM = CNSOM;
        y.litter.CNleaf = CNleaf;
        
        % Derivatives
        y.litter.dleafCdt = dleafCdt;
        y.litter.drootCdt = drootCdt;
        y.litter.dexuCdt = dexuCdt;
        y.litter.dSOMCdt = dSOMCdt;
        y.litter.dmicCrodt = dmicCdt;
        
        y.litter.drootNdt = drootNdt;
        y.litter.dNExudt = dexuCdt / params.CNE;
        y.litter.dNSOMdt = dSOMNdt;
        y.litter.dNmicrodt = dmicNdt;
        y.litter.dInorgnNdt = dInorgNdt;
        
        % Fluxes
        y.litter.A = A;
        y.litter.rootCSupply = rootCSupply;
        y.litter.leaf_resp = leaf_resp;
        y.litter.root_exu = root_exu;
        y.litter.root_resp = root_resp;
        y.litter.exu_uptake = exu_uptake;
        y.litter.SOM_uptake = SOM_uptake;
        y.litter.micro_resp = micro_resp;
        y.litter.root_N_uptake = root_N_uptake;
        y.litter.leaf_N_uptake = leaf_N_uptake;
        y.litter.phiN = phiN;
        y.litter.Imm = Imm;
        y.litter.Nmin = Nmin;
        y.litter.Nleaching = Nleaching;
        y.litter.IN_inorg = IN_inorg_to_mic;
        y.litter.micro_turnover = micro_turnover;
        y.litter.root_turnover = root_turnover;
        y.litter.rootN_uptake_norm = root_N_uptake / potential_root_N_uptake;
        
        % Growth rates
        y.litter.root_growth_rate = root_growth_rate;
        y.litter.leaf_growth_rate = leaf_growth_rate;
        y.litter.mic_growth_rate = exu_uptake + SOM_uptake - micro_resp;
        
        % Net assimilation rates
        y.litter.Net_rootC_assimilation = Net_rootC_assimilation;
        y.litter.Net_leafC_assimilation = Net_leafC_assimilation;
        
        y.litter.root_maintenance_resp = root_maintenance_resp;
        y.litter.Ccost_of_Nacq = leaf_N_uptake / rootCSupply;
        
        Anet = A - leaf_resp;
        y.litter.Anet = Anet;
        
        % For diagnostics: still compute the old constraint value
        y.litter.new_exu_constraint = ((root_maintenance_resp + root_exu) / (params.cue_R - params.min_frac_root_growth)) / Anet;
    end

end
