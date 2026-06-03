from src.state import NodeState


class DCS:
    def __init__(self, CR: float, P_line: float, q_ext: float = 0.0):
        if CR < 1.0:
            raise ValueError("CR must be >= 1.0")
        self.CR = CR
        self.P_line = P_line
        self.q_ext = q_ext

    def P_in(self) -> float:
        if self.CR == 1.0:
            return self.P_line
        return self.P_line / self.CR

    def state(self, q_std: float) -> NodeState:
        return NodeState("dcs", self.P_in(), self.P_line, self.P_line - self.P_in(), q_std, None, None, None)

