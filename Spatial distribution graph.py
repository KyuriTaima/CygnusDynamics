# Spatial distribution graph of the molecular clouds and OB star clusters in the Cygnus region

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits                   
from astropy.visualization import simple_norm

# Import csv data
objects_data = pd.read_csv('Cygnus_Objects_Datas_Schönrich_239_8.3.csv')

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
plt.figure(figsize=(12, 8))

# Add the background image
plt.imshow(img_data, extent=extent_fits, cmap='inferno', origin='lower', norm=norm_log, aspect='auto', zorder=1)

# Boucle pour ajouter les noms des objets
for i in range(len(names)):
    couleur_texte = 'blue' if 'Group' in str(names[i]) else 'black'
    if i == 4:
        xytext = (-5, -13)
    else:
        xytext = (5, 5)
    plt.annotate(
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
plt.scatter([], [], color='blue', label="OB Stars Clusters")
plt.scatter([], [], color='black', label="Molecular Clouds")
plt.scatter(l_deg, b_deg, c=point_colors)

plt.xlabel('Galactic Longitude l (degrees)')
plt.ylabel('Galactic Latitude b (degrees)')
plt.title('Spatial distribution of the molecular clouds and the OB star clusters in the Cygnus region')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.3)

# On fixe les limites d'observation
plt.xlim(83, 71)
plt.ylim(-1, 5)

plt.show()



