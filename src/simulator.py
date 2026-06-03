from __future__ import annotations

from scipy.optimize import root
import pandas as pd

from src.compressor import DCS
from src.pipe import Pipe
from src.reservoir import Reservoir
from src.state import NodeState
from src.well import Well


class FieldSimulator:
    def __init__(self, reservoir: Reservoir, wells: list[Well], shlyf: Pipe, dcs: DCS):
        if len(wells) != 3:
            raise ValueError("simulator requires exactly three wells")
        self.reservoir = reservoir
        self.wells = wells
        self.shlyf = shlyf
        self.dcs = dcs

    def solve(self, P_res: float) -> dict[str, NodeState]:
        def residuals(x):
            q1, q2, q3, P_man = x
            qs = [q1, q2, q3]
            res = []
            for q, well in zip(qs, self.wells):
                pipe_state = well.pipe.dp(P_res, max(q, 0.0))
                P_bhp = P_man + pipe_state.dP
                res.append(q - well.q(P_res, P_bhp))
            q_all = max(q1, 0.0) + max(q2, 0.0) + max(q3, 0.0) + self.dcs.q_ext
            shlyf_state = self.shlyf.dp(P_man, q_all)
            res.append(P_man - shlyf_state.dP - self.dcs.P_in())
            return res

        x0 = [500.0, 500.0, 500.0, self.dcs.P_in() + 5.0]
        solution = root(residuals, x0, method="hybr")
        q1, q2, q3, P_man = solution.x
        qs = [float(max(0.0, q1)), float(max(0.0, q2)), float(max(0.0, q3))]
        P_man = max(self.dcs.P_in(), float(P_man))

        states: dict[str, NodeState] = {}
        for i, (q, well) in enumerate(zip(qs, self.wells), start=1):
            pipe_state = well.pipe.dp(P_res, q)
            P_bhp = P_man + pipe_state.dP
            states[f"well_{i}"] = NodeState(
                f"well_{i}",
                P_bhp,
                P_man,
                pipe_state.dP,
                q,
                pipe_state.q_res,
                pipe_state.v,
                pipe_state.rho,
            )

        q_total = float(sum(qs))
        q_shlyf = q_total + self.dcs.q_ext
        shlyf_state = self.shlyf.dp(P_man, q_shlyf)
        states["shlyf"] = NodeState("shlyf", P_man, self.dcs.P_in(), shlyf_state.dP, q_shlyf, shlyf_state.q_res, shlyf_state.v, shlyf_state.rho)
        states["dcs"] = self.dcs.state(q_shlyf)
        return states

    def run(self, N_days: int, dt: float = 1.0) -> pd.DataFrame:
        rows = []
        gp = 0.0
        for step in range(N_days):
            P_res = self.reservoir.resprops.P
            states = self.solve(P_res)
            q1 = states["well_1"].q_std
            q2 = states["well_2"].q_std
            q3 = states["well_3"].q_std
            q_total = q1 + q2 + q3
            gp += q_total * dt / 1000.0
            rows.append(
                {
                    "t": step,
                    "P_res": P_res,
                    "P_man": states["shlyf"].P_in,
                    "q1": q1,
                    "q2": q2,
                    "q3": q3,
                    "q_total": q_total,
                    "Gp": gp,
                }
            )
            self.reservoir.resprops.P = self.reservoir.p2(q_total, dt)
            if step % 10 == 0 or step == N_days - 1:
                print(f"step {step + 1}/{N_days}: P_res={self.reservoir.resprops.P:.2f} atm, q_total={q_total:.1f} std.m3/day")
        return pd.DataFrame(rows)
