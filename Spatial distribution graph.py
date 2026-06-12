# Spatial distribution graph of the molecular clouds and OB star clusters in the Cygnus region

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits                   
from astropy.visualization import simple_norm
from matplotlib.patches import Ellipse # <-- NOUVEL IMPORT NÉCESSAIRE

# Import csv data
objects_data = pd.read_csv('Cygnus_Objects_Datas_Uncertainties_v3.csv')

names = objects_data['Object_Name']
distances = objects_data['Distance_pc']
l_deg = objects_data['l_deg']
b_deg = objects_data['b_deg']

# FITS image import
hdu = fits.open('Cygnus_Wise12_Grand angle_v2.fits')[0]
img_data = hdu.data

# We apply a logarithmic scale to be able to see the molecular clouds
norm_log = simple_norm(img_data, stretch='log', percent=99.5)

# Calculate the right framing (fits image is 12x12 degrees)
l_gauche = 77 + 6.0
l_droit  = 77 - 6.0
b_bas    = 2.0 - 6.0
b_haut   = 2.0 + 6.0
extent_fits = [l_gauche, l_droit, b_bas, b_haut]

# Color assignation
point_colors = ['blue' if 'Group' in str(nom) else 'black' for nom in names]

# Graph creation
fig, ax = plt.subplots(figsize=(12, 8)) # <-- Remplacé plt.figure par plt.subplots pour manipuler 'ax'

# Add the background image
ax.imshow(img_data, extent=extent_fits, cmap='inferno', origin='lower', norm=norm_log, aspect='auto', zorder=1)

# Boucle pour ajouter les noms des objets
for i in range(len(names)):
    couleur_texte = 'blue' if 'Group' in str(names[i]) else 'black'
    if i == 4:
        xytext = (-5, -13)
    else:
        xytext = (5, 5)
    ax.annotate(
        names[i], 
        (l_deg[i], b_deg[i]),          
        textcoords="offset points",    
        xytext=xytext,                 
        ha='left',                     
        fontsize=10,
        fontweight='bold',
        color=couleur_texte,
        zorder=4
    )

ax.scatter([], [], color='blue', label="OB Stars Clusters")
ax.scatter([], [], color='black', label="Molecular Clouds")
ax.scatter(l_deg, b_deg, c=point_colors, zorder=3)


# ==============================================================================
# --- ANNOTATIONS SPATIALES ---
# ==============================================================================

# 1. Ellipse Cygnus-X (autour de W75N, DR21, IRAS...)
# Tu pourras ajuster 'xy' (le centre), 'width', 'height' et 'angle' selon le rendu de ton image FITS
ellipse_cygx = Ellipse(xy=(80.8, 0.7), width=3.5, height=2.0, angle=-10, 
                       edgecolor='white', facecolor='none', linewidth=2, linestyle='--', zorder=2)
ax.add_patch(ellipse_cygx)
ax.text(80.8, 1.9, 'Cygnus-X Complex', color='white', fontsize=12, fontweight='bold', ha='center', zorder=4)

# 2. Ellipse Groupes isolés (A, C, D, F)
# Centre approximatif autour de l=74.5, b=1.3
ellipse_out = Ellipse(xy=(74.5, 1.3), width=4.5, height=1.5, angle=20, 
                      edgecolor='cyan', facecolor='none', linewidth=2, linestyle=':', zorder=2)
ax.add_patch(ellipse_out)
ax.text(74.5, 2.3, 'Outside of the cloud', color='cyan', fontsize=12, fontweight='bold', ha='center', zorder=4)

# 3. Échelle Spatiale (100 pc à 1.5 kpc)
# Calcul : 100 pc / 1500 pc = 0.0666 rad = 3.82 degrés
scale_length_deg = (100 / 1500) * (180 / np.pi)

# Position de la barre d'échelle (en bas à gauche du cadre, attention l'axe X est inversé)
l_scale_start = 76.0 
b_scale = -0.5

# Tracé de la ligne principale et des petits "bouchons" aux extrémités
ax.plot([l_scale_start, l_scale_start - scale_length_deg], [b_scale, b_scale], color='white', linewidth=3, zorder=4)
ax.plot([l_scale_start, l_scale_start], [b_scale - 0.1, b_scale + 0.1], color='white', linewidth=2, zorder=4)
ax.plot([l_scale_start - scale_length_deg, l_scale_start - scale_length_deg], [b_scale - 0.1, b_scale + 0.1], color='white', linewidth=2, zorder=4)

# Texte de la barre d'échelle
ax.text(l_scale_start - (scale_length_deg / 2), b_scale + 0.15, '100 pc (at 1.5 kpc)', 
        color='white', fontsize=10, fontweight='bold', ha='center', zorder=4)

# ==============================================================================


ax.set_xlabel('Galactic Longitude l (degrees)')
ax.set_ylabel('Galactic Latitude b (degrees)')
ax.set_title('Spatial distribution of the molecular clouds and the OB star clusters in the Cygnus region')
ax.legend(loc='upper right')
ax.grid(True, linestyle='--', alpha=0.3)

# On fixe les limites d'observation
ax.set_xlim(83, 71)
ax.set_ylim(-1, 5)

plt.show()