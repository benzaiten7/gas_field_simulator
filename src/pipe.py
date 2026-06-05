from __future__ import annotations

import math

from src.fluid import Fluid
from src.state import NodeState


class Pipe:
    # рассчитывает перепад давления при движении газа по трубе с учётом потерь на трение и гидростатитческого столба
    def __init__(self, L: float, D: float, roughness: float, fluid: Fluid, vertical_depth: float = 0.0):
        self.L = L  # длина трубы [м]
        self.D = D  # внутренний диаметр трубы [м]
        self.roughness = roughness  # шероховатость стенок [м]
        self.fluid = fluid   
        self.vertical_depth = vertical_depth  # перепад высот (глубина скважины) [м]

    def velocity(self, P: float, q: float) -> float:
        # скорость газа в трубе [м/с]
        # формула: v = 4 * Q_пласт / (π * D² * 86400),где 86400 перевод суток в секунды
        return 4.0 * self.fluid.q_res(P, q) / (math.pi * self.D * self.D * 86400.0)

    def reynolds(self, P: float, q: float) -> float:
        # число Рейнольдса (характеризует режим теения)
        # формула Re = ρ * v * D / μ
        # Re < 2300 — ламинарный режим
        # Re > 2300 — турбулентный режим
        v = self.velocity(P, q)
        mu = self.fluid.mu_pa_s(P)
        if mu <= 0.0:
            return 0.0
        return self.fluid.ro(P) * v * self.D / mu

    def friction_factor(self, Re: float) -> float:
        # для ламинарного режима (Re < 2300): λ = 64 / Re
        # для турбулентного режима: решаем уравнение Колбрука-Уайта
        # 1/√λ = -2 * log10(δ/(3.7*D) + 2.51/(Re*√λ))
        # уравнение решается методом итераций, т.е  λ подбирается за 50 попыток
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
        # расчёт перепада давления в трубе
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

