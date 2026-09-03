"""Same greenhouse OCP; light and outside temperature from forcing.csv.

CSV columns: t, light, T_out. Linear interpolants of those columns
enter f as functions of t, identical in structure to the in-code
table in greenhouse.py.
"""

import os
import sys as _sys

_sys.path.insert(0, os.path.dirname(__file__))
import greenhouse as gh

FORCING_CSV = os.path.join(os.path.dirname(__file__), "forcing.csv")

_I_FUN, _TOUT_FUN = gh.interpolants_from_file(FORCING_CSV, tag="file")
greenhouse_model = gh.make_model(_I_FUN, _TOUT_FUN)


def light(t):
    return _I_FUN(t)


def outside_temp(t):
    return _TOUT_FUN(t)


def build_ocp():
    return gh.build_ocp(model=greenhouse_model, name="greenhouse_file")


def solve(control_intervals=100, polynomial_degree=5):
    return gh.solve(
        control_intervals=control_intervals,
        polynomial_degree=polynomial_degree,
        model=greenhouse_model,
        name="greenhouse_file",
    )


if __name__ == "__main__":
    sol, sys = solve()
    print(sol)
    print("J =", sol.objective)
    print("forcing", FORCING_CSV)
    print("crop(T) =", sol.x[-1, 0], "temp(T) =", sol.x[-1, 1])
