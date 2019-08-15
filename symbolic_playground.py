import matplotlib.pyplot as plt
import numpy as np

from function_library import x, exp, log, cos

f = cos(log(exp(x) + 1))**2

x_range = np.linspace(-5, 5, 200)
N0 = 1
N = 10
for n in range(N0, N):
    f_exp = f.get_expansion(center=0.5, cut_off=n)
    plt.plot(x_range, f_exp(x_range))
plt.plot(x_range, f(x_range), c="k", ls="--", lw=2)
f_max = max(f(x_range))
f_min = min(f(x_range))
plt.ylim(f_min-0.1, f_max+0.1)
plt.show()

