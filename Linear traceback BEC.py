import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Import csv data
objects_data = pd.read_csv('Cygnus_Objects_Datas_Uncertainties.csv')

# On filtre le DataFrame pour ne garder que les Groupes B, E et C
target_names = ["Group B", "Group E", "Group C"]
supergroupBEC = objects_data[objects_data['Object_Name'].isin(target_names)].reset_index(drop=True)

# Extraction des données du repère LOCAL (Héliocentrique)
names_BEC = supergroupBEC['Object_Name'].values
X_local_pc = supergroupBEC['X_helio_pc'].values
Y_local_pc = supergroupBEC['Y_helio_pc'].values

# Extraction des vitesses PARTICULIÈRES (Locales)
U_pec_kms = supergroupBEC['U_pec_kms'].values
V_pec_kms = supergroupBEC['V_pec_kms'].values

# Time steps in million years
time_steps = [2, 5, 10, 15, 20]
# Facteur de conversion : 1 km/s = 1.0227 pc / Myr
conversion_factor = 1.0227 

# Calculate past positions
past_positions = []
for t_myr in time_steps:
    # On soustrait la distance parcourue : Pos_passée = Pos_actuelle - (Vitesse * Temps)
    X_past = X_local_pc - (U_pec_kms * conversion_factor * t_myr)
    Y_past = Y_local_pc - (V_pec_kms * conversion_factor * t_myr)
    past_positions.append((X_past, Y_past))

# Plotting
plt.figure(figsize=(10, 8))

# On trace la position actuelle d'abord
plt.scatter(X_local_pc, Y_local_pc, color='black', marker='*', s=200, label='Position Actuelle (0 Myr)', zorder=5)

for j in range(len(names_BEC)):
    plt.text(X_local_pc[j] + 10, Y_local_pc[j] + 10, names_BEC[j], fontsize=12, color='black', weight='bold')

# On trace les positions passées avec un dégradé de couleur
colors = ['#ff9999', '#ff6666', '#ff3333', '#cc0000', '#800000']
for i, t in enumerate(time_steps):
    X_past, Y_past = past_positions[i]
    plt.scatter(X_past, Y_past, label=f'-{t} Myr', color=colors[i], alpha=0.8, s=100)
    
    # Ajout de lignes pour relier visuellement la trajectoire
    if i == 0:
        for j in range(len(names_BEC)):
            plt.plot([X_local_pc[j], X_past[j]], [Y_local_pc[j], Y_past[j]], color='gray', linestyle='--', linewidth=1, alpha=0.5)
    else:
        prev_X, prev_Y = past_positions[i-1]
        for j in range(len(names_BEC)):
            plt.plot([prev_X[j], X_past[j]], [prev_Y[j], Y_past[j]], color='gray', linestyle='--', linewidth=1, alpha=0.5)

# L'axe X héliocentrique pointe vers le GC, on peut l'inverser pour que visuellement, le GC soit à gauche
plt.gca().invert_xaxis()

plt.xlabel('X local (pc) [Positif vers le Centre Galactique]')
plt.ylabel('Y local (pc) [Positif vers la rotation galactique]')
plt.title('Traceback Local Linéaire du Supergroupe BEC (Plan XY Galactique)')
plt.legend()
plt.grid(True, linestyle=':', alpha=0.7)
plt.show()