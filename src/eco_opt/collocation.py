"""Radau/Legendre orthogonal collocation coefficients.

Direct port of Yop's ``YopCollocationCoefficients`` (see
``yop-master/OptimalControl/NLP/YopCollocationCoefficients.m``), itself the
same construction as CasADi's official ``direct_collocation`` example: given
a polynomial degree ``d`` and a point selection (``'legendre'`` or
``'radau'``), builds

* ``C`` (d+1 x d+1): derivative-of-basis coefficients at each collocation
  point, used for the collocation (dynamics) equations,
* ``D`` (d+1,): basis polynomials evaluated at tau=1, used for the
  continuity/end-state equation,
* ``B`` (d+1,): integral of basis polynomials over [0,1], used for Gaussian
  quadrature of the Lagrange (running-cost) term.
"""

from __future__ import annotations

import numpy as np
import casadi as cs

__all__ = ["CollocationCoefficients"]


class CollocationCoefficients:
    def __init__(self, degree: int, points: str = "legendre"):
        self.degree = int(degree)
        self.points = points

        tau = np.array([0.0] + list(cs.collocation_points(self.degree, points)))
        d = self.degree

        C = np.zeros((d + 1, d + 1))
        D = np.zeros(d + 1)
        B = np.zeros(d + 1)
        L = np.zeros((d + 1, d + 1))  # polynomial coefficients, standard (numpy) order

        for j in range(d + 1):
            # Lagrange basis polynomial L_j(tau), built via successive convolution
            coeff = np.array([1.0])
            for r in range(d + 1):
                if r != j:
                    coeff = np.convolve(coeff, [1.0, -tau[r]])
                    coeff = coeff / (tau[j] - tau[r])
            L[j, : coeff.size] = coeff

            D[j] = np.polyval(coeff, 1.0)

            dcoeff = np.polyder(coeff)
            for r in range(d + 1):
                C[j, r] = np.polyval(dcoeff, tau[r])

            icoeff = np.polyint(coeff)
            B[j] = np.polyval(icoeff, 1.0) - np.polyval(icoeff, 0.0)

        self.tau = tau
        self.C = C
        self.D = D
        self.B = B
        self.L = L
