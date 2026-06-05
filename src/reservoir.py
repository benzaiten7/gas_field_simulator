from dataclasses import dataclass

from src.fluid import Fluid


@dataclass
class ResProps:
    P: float  # текущее пластовое давление [атм]
    V: float  # объём пласта [м³]
    T: float  # температура пласта [К]


class Reservoir:
    def __init__(self, resprops: ResProps, fluid: Fluid):
        self.resprops = resprops
        self.fluid = fluid

    def p2(self, q_total: float, dt: float = 1.0) -> float:
        # рассчитывает новое пластовое давление после отбора газа
        # аргументы: суммарный дебит всех скважин [ст.м³/сут] 
        # и шаг по времени [сут] (по умолчанию 1 сутки)
        p = self.resprops.P  
        rho = self.fluid.ro(p)
        rho_std = self.fluid.ro(self.fluid.P_STD_PA / self.fluid.ATM_TO_PA)
        z = self.fluid.z(p)
        dP = z * rho_std / rho * q_total / self.resprops.V * dt
        return max(0.0, p - dP)

