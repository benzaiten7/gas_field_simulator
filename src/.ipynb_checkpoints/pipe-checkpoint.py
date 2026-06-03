from __future__ import annotations

import math

from src.fluid import Fluid
from src.state import NodeState


class Pipe:
    def __init__(self, L: float, D: float, roughness: float, fluid: Fluid, vertical_depth: float = 0.0):
        self.L = L
        self.D = D
        self.roughness = roughness
        self.fluid = fluid
        self.vertical_depth = vertical_depth

    def velocity(self, P: float, q: float) -> float:
        return 4.0 * self.fluid.q_res(P, q) / (math.pi * self.D * self.D * 86400.0)

    def reynolds(self, P: float, q: float) -> float:
        v = self.velocity(P, q)
        mu = self.fluid.mu_pa_s(P)
        if mu <= 0.0:
            return 0.0
        return self.fluid.ro(P) * v * self.D / mu

    def friction_factor(self, Re: float) -> float:
        if Re <= 0.0:
            return 0.0
        if Re < 2300.0:
            return 64.0 / Re
        lam = 0.02
        for _ in range(50):
            rhs = -2.0 * math.log10(self.roughness / (3.7 * self.D) + 2.51 / (Re * math.sqrt(lam)))
            new_lam = 1.0 / (rhs * rhs)
            if abs(new_lam - lam) < 1e-8:
                return new_lam
            lam = new_lam
        return lam

    def dp(self, P: float, q: float) -> NodeState:
        q_abs = max(abs(q), 0.0)
        rho = self.fluid.ro(P)
        q_res = self.fluid.q_res(P, q_abs)
        v = self.velocity(P, q_abs)
        Re = self.reynolds(P, q_abs)
        lam = self.friction_factor(Re)
        friction = lam * self.L / self.D * rho * v * v / 2.0
        hydro = rho * 9.81 * self.vertical_depth
        dP = (friction + hydro) / self.fluid.ATM_TO_PA
        return NodeState("pipe", P, P - dP, dP, q, q_res, v, rho)

