from __future__ import annotations

import math

from src.fluid import Fluid
from src.pipe import Pipe


class Well:
    BETA = 0.00852702  # коэффициент пересчета единиц для измерения для формулы Дарси

    def __init__(self, fluid: Fluid, k: float, h: float, re: float, rw: float, pipe: Pipe | None = None):
        self.fluid = fluid
        self.k = k
        self.h = h
        self.re = re
        self.rw = rw
        self.pipe = pipe

    def productivity(self, P_res: float) -> float:
        mu = self.fluid.mu(P_res)
        return self.BETA * self.k * self.h / (mu * math.log(self.re / self.rw))

    def q(self, P_res: float, P_bhp: float) -> float:
        return max(0.0, self.productivity(P_res) * (P_res - P_bhp))

