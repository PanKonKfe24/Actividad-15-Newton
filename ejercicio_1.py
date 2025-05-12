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
    y_interp = np.zeros_like(x)
    for i in range(len(x)):
        term = coef[0]
        product = 1
        for j in range(1, n):
            product *= (x[i] - x_data[j-1])
            term += coef[j] * product
        y_interp[i] = term
    return y_interp

# Datos del experimento
F = np.array([50, 100, 150, 200])           # Fuerza en Newtons
epsilon = np.array([0.12, 0.35, 0.65, 1.05])  # Deformación en mm

# Estimar la deformación para 125 N
F_estimar = 125
deformacion_125 = newton_interpolation(F, epsilon, np.array([F_estimar]))[0]
print(f"La deformación estimada para una carga de {F_estimar} N es: {deformacion_125:.4f} mm")

# Graficar la interpolación
F_vals = np.linspace(min(F), max(F), 200)
epsilon_interp = newton_interpolation(F, epsilon, F_vals)

plt.figure(figsize=(8, 6))
plt.plot(F, epsilon, 'ro', label='Datos originales')
plt.plot(F_vals, epsilon_interp, 'b-', label='Interpolación de Newton')
plt.plot(F_estimar, deformacion_125, 'gs', label=f'Estimación (F={F_estimar} N)')
plt.xlabel('Carga aplicada F (N)')
plt.ylabel('Deformación ε (mm)')
plt.title('Interpolación de Newton - Deformación del material')
plt.legend()
plt.grid(True)
plt.savefig("deformacion_newton_interpolacion.png")
plt.show()
