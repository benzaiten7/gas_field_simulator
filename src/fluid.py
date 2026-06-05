from __future__ import annotations

import csv
from pathlib import Path

from src.interpolator import LinearInterpolator


class Fluid:
    R = 8.314             # универсальная газовая постоянная [Дж/(моль·К)]
    P_STD_PA = 101325.0   # стандартное давление [Па]
    T_STD = 293.15        # стандартная температура [К]
    ATM_TO_PA = 101325.0  # перевод атмосфер в Паскали

    def __init__(self, M: float, rho_c: float, xa: float, xy: float, T: float):
        self.M = M          # молярная масса газа [кг/кмоль]
        self.rho_c = rho_c  # плотность газа в стандартных условиях [кг/м³]
        self.xa = xa        # мольная доля азота
        self.xy = xy        # мольная доля диоксида углерода
        self.T = T          # температура пласта [К]
        self.x_eq = 1.0 - xa - xy # мольная доля эквивалентного углеводорода
    
        # фактор сжимаемости при стандартных условиях (формула 36 из ГОСТ)
        self.z_c = 1.0 - (0.0741 * rho_c - 0.006 - 0.063 * xa - 0.0575 * xy) ** 2
        
        # молярная масса эквивалентного углеводорода (формула 35)
        if self.x_eq > 0:
            self.M3 = (24.05525 * self.z_c * rho_c - 28.0135 * xa - 44.01 * xy) / self.x_eq
        else:
            self.M3 = 16.04  # значение по умолчанию (метан)
        
        # параметр H (формула 34)
        self.H = 128.64 + 47.479 * self.M3
        
        # загрузка таблицы вязкости
        self._mu_interpolator = self._load_mu()

    def _load_mu(self) -> LinearInterpolator:
        path = Path(__file__).resolve().parents[1] / "interp_data.csv"
        pressures = []
        viscosities = []
        
        with path.open("r", encoding="utf-8", newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                pressures.append(float(row["P_atm"]))
                viscosities.append(float(row["mu_cP"]))
        
        return LinearInterpolator(pressures, viscosities)


    def _B1(self, T: float) -> float:
        # коэффициент B1 для метана (формула 23)
        return (-0.425468 +
                2.865e-3 * T -
                4.62073e-6 * T**2 +
                (8.77118e-4 - 5.56281e-6 * T + 8.81514e-9 * T**2) * self.H +
                (-8.24747e-7 + 4.31436e-9 * T - 6.08319e-12 * T**2) * self.H**2)

    def _B2(self, T: float) -> float:
        # коэффициент B2 для азота (формула 24)
        return -0.1446 + 0.74091e-3 * T - 0.91195e-6 * T**2

    def _B23(self, T: float) -> float:
        # коэффициент B23 (формула 25)
        return -0.339693 + 0.161176e-2 * T - 0.204429e-5 * T**2

    def _B3(self, T: float) -> float:
        # коэффициент B3 для CO₂ (формула 26)
        return -0.86834 + 0.40376e-2 * T - 0.51657e-5 * T**2

    def _C1(self, T: float) -> float:
        # коэффициент C1 для метана (формула 27)
        return (-0.302488 +
                1.95861e-3 * T -
                3.16302e-6 * T**2 +
                (6.46422e-4 - 4.22876e-6 * T + 6.88157e-9 * T**2) * self.H +
                (-3.32805e-7 + 2.2316e-9 * T - 3.67713e-12 * T**2) * self.H**2)

    def _C2(self, T: float) -> float:
        # коэффициент C2 для азота (формула 28)
        return 7.8498e-3 - 3.9895e-5 * T + 6.1187e-8 * T**2

    def _C3(self, T: float) -> float:
        # коэффициент C3 для CO₂ (формула 29)
        return 2.0513e-3 + 3.4888e-5 * T - 8.3703e-8 * T**2

    def _C223(self, T: float) -> float:
        # коэффициент C223 (формула 30)
        return 5.52066e-3 - 1.68609e-5 * T + 1.57169e-8 * T**2

    def _C233(self, T: float) -> float:
        # коэффициент C233 (формула 31)
        return 3.58783e-3 + 8.06674e-6 * T - 3.25798e-8 * T**2

    def _B_star(self, T: float) -> float:
        # параметр B* (формула 32)
        return 0.72 + 1.875e-5 * (320.0 - T) ** 2

    def _C_star(self, T: float) -> float:
        # параметр C* (формула 33)
        return 0.92 + 0.0013 * (T - 270.0)


    def z(self, P: float) -> float:
        # фактор сжимаемости газа по методике GERG-91 mod. (ГОСТ 30319.2-96)
        # P: давление [атм]
        # T: температура [К], которая берётся из self.T
        T = self.T
        # Перевод давления в МПа (1 атм = 0.101325 МПа)
        P_MPa = P * 0.101325

        B1 = self._B1(T)
        B2 = self._B2(T)
        B23 = self._B23(T)
        B3 = self._B3(T)
        C1 = self._C1(T)
        C2 = self._C2(T)
        C223 = self._C223(T)
        C233 = self._C233(T)
        C3 = self._C3(T)
        B_star = self._B_star(T)
        C_star = self._C_star(T)

        # константы бинарного взаимодействия
        Y12 = 0.92
        Y13 = 0.92
        Y123 = 1.1

        # мольные доли
        x_eq = self.x_eq      # доля метана + тяжёлые
        x_N2 = self.xa        # доля азота
        x_CO2 = self.xy       # доля углекислого газа

        # расчёт Bm (формула 20 из ГОСТ)
        term1 = x_eq**2 * B1
        term2 = x_eq * x_N2 * B_star * (B1 + B2)
        term3 = -1.73 * x_eq * x_CO2 * math.sqrt(B1 * B3)
        term4 = x_N2**2 * B2
        term5 = 2.0 * x_eq * x_CO2 * B23
        term6 = x_CO2**2 * B3
        Bm = term1 + term2 + term3 + term4 + term5 + term6

        # расчёт Cm (формула 21 из ГОСТ)
        C1_2 = C1**2
        C2_2 = C2**2
        C3_2 = C3**2
        CA112 = C1_2 * C2
        CA113 = C1_2 * C3
        CA122 = C1 * C2_2
        CA123 = C1 * C2 * C3
        CA133 = C1 * C3_2

        cbrt112 = CA112 ** (1/3)
        cbrt113 = CA113 ** (1/3)
        cbrt122 = CA122 ** (1/3)
        cbrt123 = CA123 ** (1/3)
        cbrt133 = CA133 ** (1/3)

        Cm = (x_eq**3 * C1 +
              3.0 * x_eq**2 * x_N2 * cbrt112 * C_star +
              3.0 * x_eq**2 * x_CO2 * cbrt113 * Y13 +
              3.0 * x_eq * x_N2**2 * cbrt122 * C_star +
              6.0 * x_eq * x_N2 * x_CO2 * cbrt123 * Y123 +
              3.0 * x_eq * x_CO2**2 * cbrt133 * Y13 +
              x_N2**3 * C2 +
              3.0 * x_N2**2 * x_CO2 * C223 +
              3.0 * x_N2 * x_CO2**2 * C233 +
              x_CO2**3 * C3)

        # параметр b (формула 43 из ГОСТ)
        b = 1000.0 * P_MPa / (2.7715 * T)

        B0 = b * Bm
        C0 = b**2 * Cm

        A1 = 1.0 + B0
        A0 = 1.0 + 1.5 * (B0 + C0)

        # решение кубического уравнения (формулы 37-40)
        D = A0**2 - A1**3
        
        # Если D отрицательный, то возврта упрощенного значение
        if D < 0:
            return 1.0

        sqrtD = math.sqrt(D)
        A2 = (A0 - sqrtD) ** (1/3)
        z = (1.0 + A2 + A1 / A2) / 3.0

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
