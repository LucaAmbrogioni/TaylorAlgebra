import matplotlib.pyplot as plt
import numpy as np

from function_library import x, exp

f = x*exp(-x**2/(2*1))

x_range = np.linspace(-1, 1, 200)
N = 20
for n in range(2, N):
    f_exp = f.get_expansion(center=0, cut_off=n)
    plt.plot(x_range, f_exp(x_range))
plt.plot(x_range, f.evaluate(x_range), c="k", ls="--", lw=2)
plt.show()

