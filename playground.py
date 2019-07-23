import numpy as np
import matplotlib.pyplot as plt
from scipy.special import factorial

from tayloralgebra import TaylorFunction

## TODO: returns float for constant functions, to be fixed

# Single function
x_range = np.linspace(-1, 1, 100)
for K in range(2, 10):
    f = TaylorFunction([1/factorial(k) for k in range(K)])
    plt.plot(x_range, f(x_range), label=K)
plt.plot(x_range, np.exp(x_range), color="k", ls="--")
plt.show()

# Product
x_range = np.linspace(-1, 1, 100)
for K in range(2, 10):
    f = TaylorFunction([1 / factorial(k) for k in range(K)])
    g = TaylorFunction([1, 1, 3])
    h = f * g
    plt.plot(x_range, h(x_range), label=K)
plt.plot(x_range, (1 + x_range + 3 * x_range ** 2) * np.exp(x_range), color="k", ls="--")
plt.show()

# Composition #TODO: Error here (the error is probably in the composition part) (cut off??)
x_range = np.linspace(-1, 1, 100)
for K in range(2, 10):
    f = TaylorFunction([1/factorial(k) for k in range(K)], cut_off=K)
    g = TaylorFunction([1, 1, 0.5], cut_off=K)
    h = f(g)
    plt.plot(x_range, h(x_range), label=K)
plt.plot(x_range, np.exp(1 + x_range + 0.5*x_range**2), color="k", ls="--")
plt.show()

# Division #TODO: Error here
x_range = np.linspace(-1, 1, 100)
for K in range(2, 10):
    f = TaylorFunction([1/factorial(k) for k in range(K)], cut_off=K)
    g = TaylorFunction([1,1,0.5], cut_off=K)
    h = f/g
    plt.plot(x_range, h(x_range), label=K)
plt.plot(x_range, np.exp(x_range)/(1 + x_range + 0.5*x_range**2), color="k", ls="--")
plt.show()

# Everything together #TODO: Error here
x_range = np.linspace(-0.5, 0.25, 100)
for K in range(2, 10):
    f = TaylorFunction([1/factorial(k) for k in range(K)], cut_off=K)
    g = TaylorFunction([1, 1, 0.5], cut_off=K)
    h = (f/g).compose(g)
    plt.plot(x_range, h(x_range), label=K)
plt.plot(x_range, np.exp((1 + x_range + 0.5*x_range**2))/(1 + (1 + x_range + 0.5*x_range**2) + 0.5*(1 + x_range + 0.5*x_range**2)**2), color="k", ls="--")
plt.show()