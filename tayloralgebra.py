class TaylorFunction(object):

    def __init__(self, coefficients, cut_off=10):
        assert isinstance(coefficients, list) and all([isinstance(a, (float, int)) for a in coefficients])
        self.coefficients = [float(a) for a in coefficients[:cut_off]]
        self.cut_off = cut_off

    @property
    def degree(self):
        return len(self) - 1

    def __str__(self):
       return self.coefficients.__str__()

    def __call__(self, x):
        if isinstance(x, TaylorFunction):
            return self.compose(x)
        else:
            return sum([a * x ** k for k, a in enumerate(self.coefficients)])

    def __len__(self):
        return len(self.coefficients)

    def __iter__(self):
        return iter(self.coefficients)

    def __getitem__(self, index):
        return self.coefficients[index]

    def __add__(self, other):
        if isinstance(other, TaylorFunction):
            cut_off = self._get_cut_off(other)
            longer, shorter = self._sort(other)
            return TaylorFunction([a + shorter[k] if k < len(shorter) else a
                                   for k, a in enumerate(longer)], cut_off=cut_off)
        elif isinstance(other, (int, float)):
            return TaylorFunction([a + other if k == 0 else a
                                   for k, a in enumerate(self)], cut_off=self.cut_off)
        else:
            raise ValueError

    def __neg__(self):
        return -1 * self

    def __sub__(self, other):
        return self + -other

    def __radd__(self, other):
        return self + other

    def __mul__(self, other):
        if isinstance(other, (float, int)):
            cut_off = self.cut_off
            new_coefficients = [other * a for a in self.coefficients]
        elif isinstance(other, (TaylorFunction)):
            new_coefficients = []
            longer, shorter = self._sort(other)
            cut_off = self._get_cut_off(other)
            for n in range(min(len(self) + len(other) - 1, cut_off)):
                nth_coeff = sum([longer[k] * shorter[n - k] for k in range(len(longer))
                                 if len(shorter) > (n - k) >= 0])
                new_coefficients.append(nth_coeff)
        else:
            raise ValueError
        return TaylorFunction(new_coefficients, cut_off=cut_off)

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        cut_off = self._get_cut_off(other)
        inv = MultiplicativeInverse(other[0], cut_off)
        return self * inv(other - other[0])

    def compose(self, other):
        cut_off = self._get_cut_off(other)
        power = other
        result = self[0]
        for k in range(1, min(len(self), cut_off)):
            result += self[k] * power
            power = power * other
        return result

    def _sort(self, other):
        longer = self if len(self) > len(other) else other
        shorter = other if self is longer else self
        return longer, shorter

    def _get_cut_off(self, other):
        return min(self.cut_off, other.cut_off)


class MultiplicativeInverse(TaylorFunction):

    def __init__(self, a, cut_off):
        if a == 0:
            raise ZeroDivisionError
        coefficients = [(-a) ** (-n) / a for n in range(cut_off)]
        super().__init__(coefficients, cut_off)