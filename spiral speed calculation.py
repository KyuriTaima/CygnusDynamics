import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress

theta_0 = 234.8  # km/s

x_W75N = -376
y_W75N = -26
x_E = -1054
y_E = 11
x_D = -2038
y_D = 128

t_1 = 0
t_2_myr = 4.5
t_3_myr = 9.8
# Calcul des vitesses moyennes
v_W75N_E = np.sqrt((x_E - x_W75N)**2 + (y_E - y_W75N)**2) / t_2_myr
v_E_D = np.sqrt((x_D - x_E)**2 + (y_D - y_E)**2) / (t_3_myr - t_2_myr)
print(f"Vitesse moyenne de W75N à E: {v_W75N_E:.2f} pc/Myr")
print(f"Vitesse moyenne de E à D: {v_E_D:.2f} pc/Myr")
v_W75N_E_kms = v_W75N_E * 0.9778
v_E_D_kms = v_E_D * 0.9778
print(f"Vitesse moyenne de W75N à E: {v_W75N_E_kms:.2f} km/s")
print(f"Vitesse moyenne de E à D: {v_E_D_kms:.2f} km/s")

# Droite de régression linéaire pour estimer la vitesse de déplacement absolue
time_points = np.array([t_1, t_2_myr, t_3_myr])
distance_points = np.array([0, np.sqrt((x_E - x_W75N)**2 + (y_E - y_W75N)**2), 
                                np.sqrt((x_D - x_W75N)**2 + (y_D - y_W75N)**2)])
slope, intercept, r_value, p_value, std_err = linregress(time_points, distance_points)
v_regression_kms = slope * 0.9778
print(f"Vitesse de déplacement absolue estimée par régression linéaire: {v_regression_kms:.2f} km/s")
v_regression_kms_corrected = theta_0 - v_regression_kms
print(f"Vitesse de déplacement absolue corrigée par régression linéaire: {v_regression_kms_corrected:.2f} km/s")

# Correction de la vitesse de rotation galactique afin d'obtenir la vitesse de déplacement absolue
v_W75N_E_corrected = theta_0 - v_W75N_E_kms
v_E_D_corrected = theta_0 - v_E_D_kms
print(f"Vitesse corrigée de W75N à E: {v_W75N_E_corrected:.2f} km/s")
print(f"Vitesse corrigée de E à D: {v_E_D_corrected:.2f} km/s")

# Estimation de la vitesse de déplacement absolue moyenne
v_average_kms = (v_W75N_E_corrected + v_E_D_corrected) / 2
print(f"Vitesse de déplacement absolue moyenne: {v_average_kms:.2f} km/s")
