import matplotlib.pyplot as plt
import numpy as np

from function_library import x, exp

f = exp(-x**4)

x_range = np.linspace(0, 1, 200)
N0 = 1
N = 8
for n in range(N0, N):
    f_exp = f.get_expansion(center=0.5, cut_off=n)
    plt.plot(x_range, f_exp(x_range))
plt.plot(x_range, f.evaluate(x_range), c="k", ls="--", lw=2)
#plt.plot(x_range, np.exp(-x_range**4/(2*1)), c="k", ls="--", lw=2)
plt.show()

