import sys
import numpy as np
import math
import scipy.integrate


def f(x):
    f = np.empty(n)
    for i in range(n):
        xm = x[i - 1] if i - 1 >= 0 else 0
        x0 = x[i]
        xp = x[i + 1] if i + 1 < n else 0
        f[i] = xp + xm - 2 * x0 + a * ((xp - x0)**2 - (x0 - xm)**2)
    return f


def step(x0, h):
    x = x0[:n]
    v = x0[n:]

    a0 = f(x)
    x += v * h + 1 / 2 * a0 * h**2
    a1 = f(x)

    v += 1 / 2 * (a0 + a1) * h

    return np.append(x, v)


n = 32
a = 1 / 4
i = np.arange(1, n + 1)
x = np.sin(i * math.pi / n)
x = np.append(x, np.zeros(n))

sc = 4
dt = math.sqrt(1 / 8) / 4
cycles = 30000
T = 2 * math.pi / math.sqrt(3)

t = 0
tp = 0
j = 0
for j in range(cycles):
    if tp > T:
        tp -= T
        sys.stdout.write("%.16e " % (j * dt))
        for k in range(1, n + 1):
            a0 = np.dot(x[:n], np.sin(i * k * math.pi / n))
            a1 = np.dot(x[n:], np.sin(i * k * math.pi / n))
            E = 1 / 2 * a1**2 + 2 * a0**2 * math.sin(math.pi * k / (2 * n))**2
            sys.stdout.write(" %.16e" % E)
        sys.stdout.write("\n")
    t += dt
    tp += dt
    x = step(x, dt)
