from dataclasses import dataclass

from src.fluid import Fluid


@dataclass
class ResProps:
    P: float
    V: float
    T: float


class Reservoir:
    def __init__(self, resprops: ResProps, fluid: Fluid):
        self.resprops = resprops
        self.fluid = fluid

    def p2(self, q_total: float, dt: float = 1.0) -> float:
        p = self.resprops.P
        rho = self.fluid.ro(p)
        rho_std = self.fluid.ro(self.fluid.P_STD_PA / self.fluid.ATM_TO_PA)
        z = self.fluid.z(p)
        dP = z * rho_std / rho * q_total / self.resprops.V * dt
        return max(0.0, p - dP)

