import operator
import math

from scipy.special import factorial

class TaylorExpansion(object):

    def __init__(self, coefficients, center=0, cut_off=10):
        assert isinstance(coefficients, list) and all([isinstance(a, (float, int)) for a in coefficients])
        self.coefficients = [float(a) for a in coefficients[:cut_off]]
        self.center = center
        self.cut_off = cut_off

    @property
    def degree(self):
        return len(self) - 1

    def __str__(self):
       return self.coefficients.__str__()

    def __call__(self, x):
        if isinstance(x, TaylorExpansion):
            return self.compose(x)
        else:
            return sum([a * (x - self.center)**k for k, a in enumerate(self.coefficients)])

    def __len__(self):
        return len(self.coefficients)

    def __iter__(self):
        return iter(self.coefficients)

    def __getitem__(self, index):
        return self.coefficients[index]

    def __add__(self, other):
        if isinstance(other, TaylorExpansion):
            assert self.center == other.center, "You cannot sum Taylor series defined at different points (centers)"
            cut_off = self._get_cut_off(other)
            longer, shorter = self._sort(other)
            return TaylorExpansion([a + shorter[k] if k < len(shorter) else a
                                    for k, a in enumerate(longer)], center=self.center, cut_off=cut_off)
        elif isinstance(other, (int, float)):
            return TaylorExpansion([a + other if k == 0 else a
                                    for k, a in enumerate(self)], center=self.center, cut_off=self.cut_off)
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
        elif isinstance(other, (TaylorExpansion)):
            assert self.center == other.center, "You cannot multiply Taylor series defined at different points (centers)"
            new_coefficients = []
            longer, shorter = self._sort(other)
            cut_off = self._get_cut_off(other)
            for n in range(min(len(self) + len(other) - 1, cut_off)):
                nth_coeff = sum([longer[k] * shorter[n - k] for k in range(len(longer))
                                 if len(shorter) > (n - k) >= 0])
                new_coefficients.append(nth_coeff)
        else:
            raise ValueError
        return TaylorExpansion(new_coefficients, center=self.center, cut_off=cut_off)

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        cut_off = self._get_cut_off(other)
        inv = MultiplicativeInverse(other[0], cut_off)
        return self * inv(other - other[0])

    def compose(self, other):
        assert self.center == other.coefficients[0], "Error during composition: The constant coefficient of the internal function should be equal to the center of the external function"
        cut_off = self._get_cut_off(other)
        power = (other - self.center)
        result = self[0]
        for k in range(1, min(len(self), cut_off)):
            result += self[k] * power
            power = power * (other - self.center)
        return result

    def _sort(self, other):
        longer = self if len(self) > len(other) else other
        shorter = other if self is longer else self
        return longer, shorter

    def _get_cut_off(self, other):
        return min(self.cut_off, other.cut_off)


#class MultiplicativeInverse(TaylorExpansion):
    #
    # def __init__(self, a, cut_off):
    #     if a == 0:
    #         raise ZeroDivisionError
    #     coefficients = [(-a) ** (-n) / a for n in range(cut_off)]
    #     super().__init__(coefficients, cut_off)


def compose_operator(left, right):
    return left(right)


class SymbolicFunction(object):

    def get_expansion(self, center, cut_off):
        pass #Abstract

    def _apply_operator(self, other, op):
        operation_list = [op, self, other]
        return CompositeFunction(operation_list)

    def __call__(self, x):
        if isinstance(x, SymbolicFunction):
            return self.compose(x)

    def compose(self, other):
        return self._apply_operator(other, compose_operator)

    def __neg__(self):
        return -1*self

    def __add__(self, other):
        return self._apply_operator(other, operator.add)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        return self._apply_operator(other, operator.sub)

    def __rsub__(self, other):
        return -1*self.__sub__(other)

    def __mul__(self, other):
        return self._apply_operator(other, operator.mul)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        return self._apply_operator(MultiplicativeInverse()(other), operator.mul)

    def __rtruediv__(self, other):
        return self.__truediv__(other)**(-1)

    def __pow__(self, other):
        return self._apply_operator(other, operator.pow)

    def __rpow__(self, other):
        raise NotImplementedError


class PrimitiveFunction(SymbolicFunction):

    def __init__(self, coeff_generator):
        self.coeff_generator = coeff_generator
        self.taylor_expansion = {}

    def get_expansion(self, center, cut_off):
        if (center, cut_off) not in self.taylor_expansion:
            coefficients = [self.coeff_generator(center, n)
                            for n in range(cut_off)]
            self.taylor_expansion.update({(center, cut_off): TaylorExpansion(coefficients, center, cut_off)})
        return self.taylor_expansion[(center, cut_off)]


class ExponentialFunction(PrimitiveFunction):

    def __init__(self):
        coeff_generator = lambda c, n: math.exp(c)/factorial(n)
        super().__init__(coeff_generator)


class LogarithmFunction(PrimitiveFunction):

    def __init__(self):
        coeff_generator = lambda c, n: (-1)**(n+1)*c**(-n)/n if n>0 else math.log(c)
        super().__init__(coeff_generator)


class SineFunction(PrimitiveFunction):

    def __init__(self):
        coeff_generator = lambda c, n: math.sin(n*math.pi/2. + c)/factorial(n)
        super().__init__(coeff_generator)


class CosineFunction(PrimitiveFunction):

    def __init__(self):
        coeff_generator = lambda c, n: math.cos(n*math.pi/2. + c)/factorial(n)
        super().__init__(coeff_generator)


class MultiplicativeInverse(PrimitiveFunction):

    def __init__(self):
        coeff_generator = lambda c, n: (-1)**n*c**(-(n+1))
        super().__init__(coeff_generator)


class PowerFunction(PrimitiveFunction):

    def __init__(self, k):
        def decreasing_factorial(k, n):
            if n < 0:
                return 0
            elif n == 0:
                return 1
            elif n == 1:
                return k
            else:
                return decreasing_factorial(k, n-1)
        coeff_generator = lambda c, n: (decreasing_factorial(k, n)/factorial(n))*c**(k-n) if n <= k else 0.
        super().__init__(coeff_generator)

class CompositeFunction(SymbolicFunction):

    def __init__(self, operation_list):
        self.operation_list = operation_list

    def get_expansion(self, center, cut_off):
        op, left, right = self.operation_list
        if op is compose_operator:
            left_center = right.get_expansion(center, cut_off)[0]
        else:
            left_center = center
        return op(left.get_expansion(left_center, cut_off) if isinstance(left, SymbolicFunction) else left,
                  right.get_expansion(center, cut_off) if isinstance(right, SymbolicFunction) else right)
