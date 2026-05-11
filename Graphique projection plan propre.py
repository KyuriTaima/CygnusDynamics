# Graphic production for analysis of molecular clouds and OB star clusters in Cygnus X

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Import csv data for molecular clouds and OB star clusters from a csv file
objects_data = pd.read_csv('Cygnus_Objects_Datas_Reid2019_236_8.15.csv')
# Extract relevant columns for clouds and clusters
names = objects_data['Object_Name']
U_pec_kms = objects_data['U_pec_kms']
V_pec_kms = objects_data['V_pec_kms']
W_pec_kms = objects_data['W_pec_kms']
radial_velocity_kms = objects_data['Radial_Velocity_kms']
v_l_kms = objects_data['v_l_kms']
v_b_kms = objects_data['v_b_kms']
distances = objects_data['Distance_pc']
l_deg = objects_data['l_deg']
b_deg = objects_data['b_deg']
mu_l = objects_data['mu_l_mas_yr']
mu_b = objects_data['mu_b_mas_yr']

# Calculation of the velocity vector in the observation plane (v_l, v_b) from the proper motions and distances
v = np.sqrt(v_l_kms**2 + v_b_kms**2)

# Color assignation based on type of object
arrow_colors = ['red' if 'Group' in str(nom) else 'cyan' for nom in names]

# Plotting the velocity vectors in the observation plane (v_l, v_b) for molecular clouds and OB star clusters
plt.figure(figsize=(10, 8))
plt.quiver(l_deg, b_deg, v_l_kms, v_b_kms, angles='xy', scale_units='xy', scale=5, color=arrow_colors, width = 0.003, headwidth = 3.5, headlength = 4)

# Boucle pour ajouter les noms des objets
for i in range(len(names)):
    plt.annotate(
        names[i], 
        (l_deg[i], b_deg[i]),          # Coordonnées du point d'origine
        textcoords="offset points",    # On indique qu'on veut décaler le texte
        xytext=(5, 5),                 # Décalage de 5 pixels à droite et 5 en haut
        ha='left',                     # Alignement horizontal
        fontsize=9,
        fontweight='bold',
        color='darkred'                # Une couleur différente pour bien distinguer le texte
    )

plt.xlabel('Galactic Longitude (degrees)')
plt.ylabel('Galactic Latitude (degrees)')
plt.title('Velocity Vectors of Molecular Clouds and OB Star Clusters in Cygnus X')
plt.grid()
plt.xlim(83, 70)
plt.ylim(-1, 5)
plt.show() 


