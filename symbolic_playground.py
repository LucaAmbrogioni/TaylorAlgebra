import matplotlib.pyplot as plt
import numpy as np

from tayloralgebra import ExponentialFunction, LogarithmFunction, SineFunction, PowerFunction

A = LogarithmFunction()

x_range = np.linspace(1, 3, 200)
for N in range(2, 10):
    f = A.get_expansion(center=2, cut_off=N)
    plt.plot(x_range, f(x_range))
plt.plot(x_range, np.log(x_range), c="k")
plt.show()

