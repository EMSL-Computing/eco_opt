"""Optimality-based plant carbon-nitrogen allocation model
(Chakrawal, Mooney & Colombi 2025/2026, "Modelling belowground plant
acclimation to low soil nitrogen", New Phytologist,
https://doi.org/10.1111/nph.71210; code DOI
https://doi.org/10.5281/zenodo.17872549).

Nested/hierarchical control: exudation is a *fraction* of net root
carbon assimilation rather than an independent absolute rate, which
avoids artificial competition between root growth and exudation.

Problem (fixed terminal time, "base" scenario)::

    maximize   int_0^T leaf_growth_rate dt
    subject to  9-state plant-soil-microbe C/N dynamics (see
                ``plant_soil_microbe_nested`` below), 2 controls:
                  u0 = rootCSupply       [gC/m2/d]  (leaf->root C allocation)
                  u1 = exu_fraction      [-]        (fraction of net root C
                                                      assimilation exuded)
                state/path bounds on C/N pools, growth rates, uptake ratio,
                and allocation fractions (see ``build_ocp``)
                x(0) = state_init (leafC=50, rootC=75, exuC=0.5, SOMC=7200,
                micC=2.8, plus stoichiometric N pools, inorgN=20)
                t_f == terminal_time (120 d default)

This is the full 9-state / 2-control ecological OCP underlying the
project's earlier (simpler) examples; unlike those, no closed-form or
external reference solution exists, so this port is validated via
(a) `Solve_Succeeded` convergence, (b) non-negativity / mass-balance
sanity of all 9 states, and (c) a regression-tracked objective value.
"""

import casadi as cs

import eco_opt
from eco_opt.system import Namespace

# ---------------------------------------------------------------------------
# Parameters ("base" scenario)
# ---------------------------------------------------------------------------
PARAMS = dict(
    a=600 / 365,          # max nitrogen productivity gC/gN.d
    b=0.001,               # decreasing factor in P_N = a - b*C_P
    KNL=75e8,               # max leaf N (effectively disables logistic term)
    Vmax_rootCsupply=1.5,   # (unused directly in dynamics; retained for parity)
    Vmax_exu_production=0.0019 * 4,  # (unused directly; retained for parity)
    tau_root=0.4 / 365,     # root turnover rate constant d^-1
    Vmax_ex=5.3 / 10,       # max. exudate uptake rate d^-1
    KM_ex=3.3,              # half-sat const. for exudate uptake gC/m2
    Vmax_SOM=0.15,          # max. SOM uptake rate d^-1
    KM_SOM=6443.0,          # half-sat const. for SOM uptake rate gC/m2
    tau_mic=5 / 365,        # microbial turnover rate constant d^-1
    Vmax_IN_to_root=3.7 / 15,  # max. root inorganic N uptake rate d^-1
    KM_IN_to_root=20.0,     # half-sat const. for root inorganic N uptake
    Vmax_IN_to_mic=15.0,    # max. inorganic N supply rate to microbial pool d^-1
    KM_IN_to_mic=1.0,       # half-sat const. for inorganic N supply to microbes gN/m2
    leaching_rate=0.0,      # inorganic N leaching rate d^-1
    e_L=0.7,                # leaf C assimilation efficiency
    cue_R=0.6,              # root carbon-use efficiency
    cue_M=0.3,              # microbial carbon-use efficiency
    CNE=10000.0,            # C:N ratio of exudates
    CNM_clim=10.0,          # C-limited C:N ratio of microbial pool
    CNleaf_opt=21.0,        # optimal C:N ratio of leaf
    CNroot_opt=50.0,        # C:N ratio of roots
    r_m=0.005,              # root maintenance respiration rate constant
    min_frac_leaf_growth=0.2,
    min_frac_root_growth=0.2,
    min_exudation=0.1,
    max_exudation=0.5,
)

# Initial state ("base" scenario)
_LEAFC0 = 50.0
_ROOTC0 = 1.5 * _LEAFC0
_CNSOM0 = 10.0
_SOMC0 = 7200.0
STATE_INIT = dict(
    leafC=_LEAFC0,
    rootC=_ROOTC0,
    exuC=0.5,
    SOMC=_SOMC0,
    micC=2.8,
    leafN=_LEAFC0 / PARAMS["CNleaf_opt"],
    rootN=_ROOTC0 / PARAMS["CNroot_opt"],
    SOMN=_SOMC0 / _CNSOM0,
    micN=2.8 / PARAMS["CNM_clim"],
    inorgN=20.0,
)


def plant_soil_microbe_nested(t, state, u, params=PARAMS):
    """9-state plant-soil-microbe C/N dynamics, nested control structure.

    State order:
        0 leafC, 1 rootC, 2 exuC, 3 SOMC, 4 micC,
        5 rootN, 6 SOMN, 7 micN, 8 inorgN
    Control order: 0 rootCSupply [gC/m2/d], 1 exu_fraction [-]
    """
    p = params
    fCL, fCR, fCE, fSOMC, fCM = state[0], state[1], state[2], state[3], state[4]
    fNR, fNSOM, fNM, fNinorg = state[5], state[6], state[7], state[8]
    fNL = fCL / p["CNleaf_opt"]

    rootCSupply = u[0]
    exu_fraction = u[1]

    CNR = fCR / fNR
    CNM = fCM / fNM
    CNSOM = fSOMC / fNSOM
    CNleaf = fCL / fNL

    A = fNL * (1 - fNL / p["KNL"]) * (p["a"] - p["b"] * fCL)
    leaf_resp = (1 - p["e_L"]) * A
    root_turnover = p["tau_root"] * fCR
    exu_uptake = (p["Vmax_ex"] * fCE * fCM) / (p["KM_ex"] + fCM + fCE)
    SOM_uptake = (p["Vmax_SOM"] * fSOMC * fCM) / (p["KM_SOM"] + fCM + fSOMC)
    micro_turnover = p["tau_mic"] * fCM
    micro_resp = (1 - p["cue_M"]) * (exu_uptake + SOM_uptake)
    potential_root_N_uptake = (p["Vmax_IN_to_root"] * fNR * fNinorg) / (p["KM_IN_to_root"] + fNinorg)

    root_resp = (1 - p["cue_R"]) * rootCSupply
    root_maintenance_resp = p["r_m"] * fCR
    leaf_N_uptake = (A - leaf_resp - rootCSupply) / p["CNleaf_opt"]

    Net_rootC_assimilation = rootCSupply - root_resp
    root_exu = exu_fraction * Net_rootC_assimilation
    root_growth_rate = Net_rootC_assimilation - root_exu - root_maintenance_resp

    phiC = (SOM_uptake / CNSOM + exu_uptake / p["CNE"]
            - p["cue_M"] * (exu_uptake + SOM_uptake) / CNM)
    IN_inorg_to_mic = p["Vmax_IN_to_mic"] * fNinorg / (p["KM_IN_to_mic"] + fNinorg)
    Nleaching = p["leaching_rate"] * fNinorg

    Nmin = cs.if_else(phiC > 0, phiC, 0)
    Imm = cs.if_else(phiC > 0, 0, cs.if_else(-phiC > IN_inorg_to_mic, IN_inorg_to_mic, -phiC))
    phiN = Nmin - Imm

    Net_leafC_assimilation = A - leaf_resp
    leaf_growth_rate = Net_leafC_assimilation - rootCSupply
    frac_leaf_growth = leaf_growth_rate / Net_leafC_assimilation
    frac_root_growth = root_growth_rate / Net_rootC_assimilation
    root_N_uptake = root_growth_rate / p["CNroot_opt"] + leaf_N_uptake + root_exu / p["CNE"]
    rootN_uptake_norm = root_N_uptake / potential_root_N_uptake

    dleafC = leaf_growth_rate
    drootC = root_growth_rate - root_turnover
    dexuC = root_exu - exu_uptake
    dSOMC = -SOM_uptake + micro_turnover + root_turnover
    dmicC = exu_uptake + SOM_uptake - micro_resp - micro_turnover
    drootN = root_N_uptake - root_exu / p["CNE"] - root_turnover / p["CNroot_opt"] - leaf_N_uptake
    dSOMN = -SOM_uptake / CNSOM + micro_turnover / CNM + root_turnover / p["CNroot_opt"]
    dmicN = exu_uptake / p["CNE"] + SOM_uptake / CNSOM - micro_turnover / CNM - phiN
    dInorgN = phiN - root_N_uptake - Nleaching

    dx = cs.vertcat(dleafC, drootC, dexuC, dSOMC, dmicC, drootN, dSOMN, dmicN, dInorgN)

    y = Namespace(
        leafC=fCL, rootC=fCR, exuC=fCE, SOMC=fSOMC, micC=fCM,
        rootN=fNR, SOMN=fNSOM, micN=fNM, inorgN=fNinorg, leafN=fNL,
        CNR=CNR, CNM=CNM, CNSOM=CNSOM, CNleaf=CNleaf,
        A=A, Anet=Net_leafC_assimilation, leaf_resp=leaf_resp,
        rootCSupply=rootCSupply, exu_fraction=exu_fraction,
        root_exu=root_exu, root_resp=root_resp,
        root_maintenance_resp=root_maintenance_resp,
        Net_rootC_assimilation=Net_rootC_assimilation,
        Net_leafC_assimilation=Net_leafC_assimilation,
        leaf_growth_rate=leaf_growth_rate, root_growth_rate=root_growth_rate,
        frac_leaf_growth=frac_leaf_growth, frac_root_growth=frac_root_growth,
        leaf_N_uptake=leaf_N_uptake, root_N_uptake=root_N_uptake,
        rootN_uptake_norm=rootN_uptake_norm,
        exu_uptake=exu_uptake, SOM_uptake=SOM_uptake,
        micro_turnover=micro_turnover, root_turnover=root_turnover,
        phiN=phiN, phiC=phiC,
    )
    return dx, y


def build_ocp(terminal_time=120.0, state_init=STATE_INIT, params=PARAMS):
    sys = eco_opt.System(
        states=9, controls=2,
        model=lambda t, x, u: plant_soil_microbe_nested(t, x, u, params),
        name="root_opt",
    )
    time = sys.t
    litter = sys.y

    ocp = eco_opt.Ocp()
    ocp.max([litter.leaf_growth_rate])
    ocp.st(
        system=sys,
        constraints=[
            (0, "<=", litter.leafC),
            (0, "<=", litter.rootC),
            (0, "<=", litter.exuC),
            (0, "<=", litter.SOMC),
            (0, "<=", litter.micC),
            (0, "<=", litter.rootN),
            (0, "<=", litter.SOMN),
            (0, "<=", litter.micN),
            (0.01, "<=", litter.inorgN),
            (0, "<=", litter.CNM, "<=", 20),
            (0, "<=", litter.leaf_growth_rate),
            (0, "<=", litter.root_exu),
            (0, "<=", litter.rootCSupply),
            (params["min_exudation"], "<=", sys.u[1], "<=", params["max_exudation"]),
            (0, "<=", litter.rootN_uptake_norm, "<=", 1),
            (params["min_frac_leaf_growth"], "<=", litter.frac_leaf_growth),
            (params["min_frac_root_growth"], "<=", litter.frac_root_growth),
            (0, "==", eco_opt.t_0(time)),
            (state_init["leafC"], "==", eco_opt.t_0(litter.leafC)),
            (state_init["rootC"], "==", eco_opt.t_0(litter.rootC)),
            (state_init["exuC"], "==", eco_opt.t_0(litter.exuC)),
            (state_init["SOMC"], "==", eco_opt.t_0(litter.SOMC)),
            (state_init["micC"], "==", eco_opt.t_0(litter.micC)),
            (state_init["rootN"], "==", eco_opt.t_0(litter.rootN)),
            (state_init["SOMN"], "==", eco_opt.t_0(litter.SOMN)),
            (state_init["micN"], "==", eco_opt.t_0(litter.micN)),
            (state_init["inorgN"], "==", eco_opt.t_0(litter.inorgN)),
            (terminal_time, "==", eco_opt.t_f(time)),
        ],
    )
    return ocp, sys


def solve(terminal_time=120.0, state_init=STATE_INIT, params=PARAMS,
          control_intervals=50, polynomial_degree=3):
    ocp, sys = build_ocp(terminal_time=terminal_time, state_init=state_init, params=params)
    # Reasonable initial guesses matching the cold-start scale
    # (states in the hundreds/thousands, controls O(0.1-1)) help ipopt
    # converge in far fewer iterations than the package default of 1.0
    # everywhere.
    init_guess = {
        sys.x[0]: state_init["leafC"],
        sys.x[1]: state_init["rootC"],
        sys.x[2]: state_init["exuC"],
        sys.x[3]: state_init["SOMC"],
        sys.x[4]: state_init["micC"],
        sys.x[5]: state_init["rootN"],
        sys.x[6]: state_init["SOMN"],
        sys.x[7]: state_init["micN"],
        sys.x[8]: state_init["inorgN"],
        sys.u[0]: 1.0,
        sys.u[1]: 0.25,
    }
    sol = ocp.solve(
        control_intervals=control_intervals,
        polynomial_degree=polynomial_degree,
        collocation_points="legendre",
        initial_guess=init_guess,
        ipopt_options={"max_iter": 2000, "tol": 1e-6},
        tf_guess=terminal_time,
    )
    return sol, sys


# Regression-tracked value: computed with eco_opt at the base scenario
# above, control_intervals=50, polynomial_degree=3, legendre collocation.
# Updated when path inequalities were also applied at interval starts.
REGRESSION_OBJECTIVE = 499.3902440480594


if __name__ == "__main__":
    sol, sys = solve()
    print(sol)
    print("objective (cumulative leaf growth):", sol.objective)
    print("leafC(0)->leafC(tf):", sol.x[0, 0], "->", sol.x[-1, 0])
    print("rootC(0)->rootC(tf):", sol.x[0, 1], "->", sol.x[-1, 1])
    print("inorgN(0)->inorgN(tf):", sol.x[0, 8], "->", sol.x[-1, 8])
