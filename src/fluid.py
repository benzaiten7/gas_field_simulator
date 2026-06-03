from __future__ import annotations

import csv
from pathlib import Path

from src.interpolator import LinearInterpolator


class Fluid:
    R = 8.314
    P_STD_PA = 101325.0
    T_STD = 293.15
    ATM_TO_PA = 101325.0

    def __init__(self, M: float, rho_c: float, xa: float, xy: float, T: float):
        self.M = M
        self.rho_c = rho_c
        self.xa = xa
        self.xy = xy
        self.T = T
        self._mu_interpolator = self._load_mu()

    def _load_mu(self) -> LinearInterpolator:
        path = Path(__file__).resolve().parents[1] / "interp_data.csv"
        ps: list[float] = []
        mus: list[float] = []
        with path.open("r", encoding="utf-8", newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                ps.append(float(row["P_atm"]))
                mus.append(float(row["mu_cP"]))
        return LinearInterpolator(ps, mus)

    def z(self, P: float) -> float:
        p = max(P, 0.0)
        composition = 1.0 + 0.08 * self.xy - 0.03 * self.xa
        z = 1.0 - 0.165 * composition * p / (p + 27.0) + 0.0000015 * p * p
        return max(0.65, min(1.25, z))

    def ro(self, P: float) -> float:
        return P * self.ATM_TO_PA * self.M / (self.z(P) * self.R * self.T)

    def bg(self, P: float) -> float:
        return self.P_STD_PA * self.z(P) * self.T / (P * self.ATM_TO_PA * self.T_STD)

    def mu(self, P: float) -> float:
        return self._mu_interpolator.predict(P)

    def mu_pa_s(self, P: float) -> float:
        return self.mu(P) / 1000.0

    def q_res(self, P: float, q_std: float) -> float:
        return q_std * self.bg(P)

