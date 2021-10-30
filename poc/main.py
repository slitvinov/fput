import sys
import numpy as np
import math
import scipy.integrate


def f(x, t):
    f = np.empty(2 * n)
    for i in range(n):
        f[i] = x[i + n]
    for i in range(n):
        xm = x[i - 1] if i - 1 >= 0 else 0
        x0 = x[i]
        xp = x[i + 1] if i + 1 < n else 0
        f[i + n] = xp + xm - 2 * x0 + a * ((xp - x0)**2 - (x0 - xm)**2)
    return f


n = 32
a = 1 / 4
i = np.arange(1, n + 1)
x0 = np.sin(i * math.pi / n)
x0 = np.append(x0, np.zeros(n))

dt = math.sqrt(1 / 8)
T = 30000 * dt
t = np.linspace(0, T, 2000)
x = scipy.integrate.odeint(f, x0, t)

for j in range(len(x)):
    for k in range(1, n + 1):
        a0 = np.dot(x[j][:n], np.sin(i * k * math.pi / n))
        a1 = np.dot(x[j][n:], np.sin(i * k * math.pi / n))
        E = 1 / 2 * a1**2 + 2 * a0**2 * math.sin(math.pi * k / (2 * n))**2
        sys.stdout.write("%.16e " % E)
    sys.stdout.write("\n")
