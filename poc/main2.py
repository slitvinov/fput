import sys
import numpy as np
import math
import scipy.integrate


def f(x):
    f = np.empty(2 * n)
    for i in range(n):
        f[i] = x[i + n]
    for i in range(n):
        xm = x[i - 1] if i - 1 >= 0 else 0
        x0 = x[i]
        xp = x[i + 1] if i + 1 < n else 0
        f[i + n] = xp + xm - 2 * x0 + a * ((xp - x0)**2 - (x0 - xm)**2)
    return f


def rk(x, h):
    k1 = f(x)
    k2 = f(x + h * k1 / 2)
    k3 = f(x + h * k2 / 2)
    k4 = f(x + h * k3)
    return x + 1 / 6 * h * (k1 + 2 * k2 + 2 * k3 + k4)


n = 32
a = 1 / 4
i = np.arange(1, n + 1)
x = np.sin(i * math.pi / n)
x = np.append(x, np.zeros(n))

j = 0
t = 0
tp = 0
w = 2 * math.sin(math.pi / (2 * (n + 1)))
T = 2 * math.pi / w
dt = T / 100
while True:
    if t == 0 or tp > T:
        tp -= T
        sys.stdout.write("%.16e " % t)
        for k in range(1, n + 1):
            a0 = np.dot(x[:n], np.sin(i * k * math.pi / n))
            a1 = np.dot(x[n:], np.sin(i * k * math.pi / n))
            E = 1 / 2 * a1**2 + 2 * a0**2 * math.sin(math.pi * k / (2 * n))**2
            sys.stdout.write(" %.16e" % E)
        sys.stdout.write("\n")
    t += dt
    tp += dt
    if t > 30000 * math.sqrt(1 / 8):
        break
    x = rk(x, dt)
