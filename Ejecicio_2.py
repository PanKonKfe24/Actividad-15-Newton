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
temperaturas = np.array([200, 250, 300, 350, 400])  # T (°C)
eficiencias = np.array([30, 35, 40, 46, 53])         # Eficiencia (%)

# Predicción para T = 275°C
T_objetivo = 275
eficiencia_predicha = newton_interpolation(temperaturas, eficiencias, np.array([T_objetivo]))[0]

print(f"Eficiencia estimada para T = {T_objetivo}°C: {eficiencia_predicha:.2f}%")

# Gráfica
x_vals = np.linspace(min(temperaturas), max(temperaturas), 200)
y_interp = newton_interpolation(temperaturas, eficiencias, x_vals)

plt.figure(figsize=(8, 6))
plt.plot(temperaturas, eficiencias, 'ro', label='Datos experimentales')
plt.plot(x_vals, y_interp, 'b-', label='Interpolación de Newton')
plt.plot(T_objetivo, eficiencia_predicha, 'gs', label=f'Predicción ({T_objetivo}°C)')
plt.xlabel('Temperatura (°C)')
plt.ylabel('Eficiencia (%)')
plt.title('Interpolación de Newton - Motor térmico')
plt.legend()
plt.grid(True)
plt.savefig("eficiencia_motor_newton.png")
plt.show()
