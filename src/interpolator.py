class LinearInterpolator:
    def __init__(self, xs: list, ys: list):
        if len(xs) != len(ys):
            raise ValueError("xs and ys must have the same length")
        if len(xs) < 2:
            raise ValueError("at least two points are required")
        for i in range(len(xs) - 1):
            if xs[i] >= xs[i + 1]:
                raise ValueError("xs must be sorted in ascending order")
        self.xs = list(xs)
        self.ys = list(ys)

    def predict(self, xp: float) -> float:
        if xp < self.xs[0] or xp > self.xs[-1]:
            raise ValueError("point is outside interpolation range")
        for i in range(len(self.xs) - 1):
            x0 = self.xs[i]
            x1 = self.xs[i + 1]
            if x0 <= xp <= x1:
                y0 = self.ys[i]
                y1 = self.ys[i + 1]
                return y0 + (y1 - y0) / (x1 - x0) * (xp - x0)
        return self.ys[-1]

