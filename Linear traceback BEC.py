import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Import csv data
objects_data = pd.read_csv('Cygnus_Objects_Datas_Uncertainties.csv')

# On filtre le DataFrame pour ne garder que les Groupes B, E et C
target_names = ["Group B", "Group E", "Group C"]
supergroupBEC = objects_data[objects_data['Object_Name'].isin(target_names)].reset_index(drop=True)

# On extrait les noms et on convertit les colonnes en tableaux NumPy
names_BEC = supergroupBEC['Object_Name'].values
X_gc_kpc = supergroupBEC['X_gc_pc'].values
Y_gc_kpc = supergroupBEC['Y_gc_pc'].values
Z_gc_kpc = supergroupBEC['Z_gc_pc'].values
vx_gc_kms = supergroupBEC['vx_gc_kms'].values
vy_gc_kms = supergroupBEC['vy_gc_kms'].values
vz_gc_kms = supergroupBEC['vz_gc_kms'].values

# Time steps in million years
time_steps = [2, 5, 10, 15, 20]
time_steps_sec = [t * 1e6 * 365.25 * 24 * 3600 for t in time_steps]

# Calculate past positions
past_positions = []
for t_sec in time_steps_sec:
    X_past = X_gc_kpc - vx_gc_kms * t_sec / (3.086e+16)
    Y_past = Y_gc_kpc - vy_gc_kms * t_sec / (3.086e+16)
    Z_past = Z_gc_kpc - vz_gc_kms * t_sec / (3.086e+16)
    past_positions.append((X_past, Y_past, Z_past))

# Plotting
plt.figure(figsize=(10, 8))
# On trace la position actuelle d'abord
plt.scatter(X_gc_kpc, Y_gc_kpc, color='black', marker='*', s=150, label='Current Position')
for j in range(len(names_BEC)):
    plt.text(X_gc_kpc[j], Y_gc_kpc[j], names_BEC[j], fontsize=10, color='black', weight='bold')

# On trace le passé
for i, t in enumerate(time_steps):
    X_past, Y_past, Z_past = past_positions[i]
    plt.scatter(X_past, Y_past, label=f'-{t} Myr', alpha=0.7)
    
plt.xlabel('X (kpc)')
plt.ylabel('Y (kpc)')
plt.title('Linear Traceback of Supergroup BEC (XY Plane)')
plt.legend()
plt.grid()
plt.show()