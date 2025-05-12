# -*- coding: utf-8 -*-
#           Autor:
#   Buitrago Rios Juan Pablo
#   
#   Fecha: 01/05/2005
#   Version: 1.01

import numpy as np
import matplotlib.pyplot as plt

def newton_divided_diff(x, y):
    n = len(x)
    coef = np.zeros([n, n])
    coef[:, 0] = y

    for j in range(1, n):
        for i in range(n - j):
            coef[i, j] = (coef[i+1, j-1] - coef[i, j-1]) / (x[i+j] - x[i])

    return coef[0, :]

def newton_interpolation(x_data, y_data, x):
    coef = newton_divided_diff(x_data, y_data)
    n = len(x_data)

    y_interp = np.zeros_like(x, dtype=float)
    for i in range(len(x)):
        term = coef[0]
        product = 1
        for j in range(1, n):
            product *= (x[i] - x_data[j-1])
            term += coef[j] * product
        y_interp[i] = term

    return y_interp

# Datos experimentales
velocidades = np.array([10, 20, 30, 40, 50, 60])   # V (m/s)
coef_arrastre = np.array([0.32, 0.30, 0.28, 0.27, 0.26, 0.25])  # Cd

# Predicción para V = 35 m/s
V_objetivo = 35
Cd_estimado = newton_interpolation(velocidades, coef_arrastre, np.array([V_objetivo]))[0]

print(f"Coeficiente de arrastre estimado para V = {V_objetivo} m/s: {Cd_estimado:.4f}")

# Gráfica
x_vals = np.linspace(min(velocidades), max(velocidades), 200)
y_interp = newton_interpolation(velocidades, coef_arrastre, x_vals)

plt.figure(figsize=(8, 6))
plt.plot(velocidades, coef_arrastre, 'ro', label='Datos reales')
plt.plot(x_vals, y_interp, 'b-', label='Interpolación de Newton')
plt.plot(V_objetivo, Cd_estimado, 'gs', label=f'Estimación ({V_objetivo} m/s)')
plt.xlabel('Velocidad del aire (m/s)')
plt.ylabel('Coeficiente de arrastre $C_d$')
plt.title('Interpolación de Newton - Coeficiente de arrastre')
plt.legend()
plt.grid(True)
plt.savefig("coef_arrastre_newton.png")
plt.show()
