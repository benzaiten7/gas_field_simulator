class LinearInterpolator:
    #Класс, который по таблмце значений (x,y) находит y для любого x внутри диапазона
    def __init__(self, xs: list, ys: list):
        #Сохраняет таблицу для интерполяции, аргумент xs это список значений x и он должен быть отсортирован по возрастанию
        if len(xs) != len(ys):
            raise ValueError("xs и ys должны быть одинаковой длины")
        if len(xs) < 2: 
            raise ValueError("нужно минимум 2 точки для интерполяции")
        for i in range(len(xs) - 1):
            if xs[i] >= xs[i + 1]:
                raise ValueError("xs должен быть отсортирован по возрастанию")
        self.xs = list(xs)
        self.ys = list(ys)

    def predict(self, xp: float) -> float:
        #Возвращает y для заданного xp методом линейной интерполяции
        #Если xp находится между двумя точками таблицы, вычисляет y по формуле прямой
        #Если xp совпадает с точкой из таблицы, возвращает соответствующий y
        if xp < self.xs[0] or xp > self.xs[-1]:
            raise ValueError("Точка вне диапазона таблицы")
        for i in range(len(self.xs) - 1):
            #Поиск между какими двумя точками находится xp
            x0 = self.xs[i]
            x1 = self.xs[i + 1]
            if x0 <= xp <= x1: #Если xp между x0 b x1 (или совпадает с одним из них)
                y0 = self.ys[i]
                y1 = self.ys[i + 1]
                return y0 + (y1 - y0) / (x1 - x0) * (xp - x0)
        return self.ys[-1] #Если ничего не нашли

