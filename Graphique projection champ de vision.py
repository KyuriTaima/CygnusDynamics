# Graphic production for analysis of molecular clouds and OB star clusters in Cygnus X

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.io import fits                   
from astropy.visualization import simple_norm

# Import csv data
objects_data = pd.read_csv('Cygnus_Objects_Datas_Reid2019_236_8.15.csv')

# Extract relevant columns
names = objects_data['Object_Name'].values
U_pec_kms = objects_data['U_pec_kms'].values
V_pec_kms = objects_data['V_pec_kms'].values
W_pec_kms = objects_data['W_pec_kms'].values
l_deg = objects_data['l_deg'].values
b_deg = objects_data['b_deg'].values

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


# Graph creation
plt.figure(figsize=(12, 8))

# Affichage du fond étoilé/nuageux en premier (zorder=1)
plt.imshow(img_data, extent=extent_fits, cmap='inferno', origin='lower', norm=norm_log, aspect='auto', zorder=1)

# Couleurs des flèches
arrow_colors = ['blue' if 'Group' in str(nom) else 'black' for nom in names]

# Superposition des vecteurs vitesse (zorder=3 pour être au-dessus du fond)
plt.quiver(l_deg, b_deg, -U_pec_kms, W_pec_kms, 
           angles='xy', scale_units='xy', scale=10, 
           color=arrow_colors, width=0.003, headwidth=3.5, headlength=4, zorder=3)

# Boucle pour ajouter les noms des objets
for i in range(len(names)):
    # On met le texte de la même couleur que la flèche pour bien lire sur le fond sombre
    couleur_texte = 'blue' if 'Group' in str(names[i]) else 'black'
    
    if i == 6: # On met un décalage différent pour le dernier objet (DR20) pour éviter qu'il ne soit superposé à DR21
        xytext = (-35, -10)
    elif i== 8:
        xytext = (0, -10)
    elif i == 4:
        xytext = (-30, 5)
    elif i == 1:
        xytext = (0, -10)
    else:
        xytext = (3, 0)
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

plt.xlabel('Galactic Longitude l (degrees)')
plt.ylabel('Galactic Latitude b (degrees)')
plt.title('Velocity Vectors (-U, W) overlaid on WISE 12 Image of Cygnus X')

plt.grid(True, linestyle='--', alpha=0.3)

# On fixe les limites d'observation
plt.xlim(83, 70)
plt.ylim(-1, 5)

plt.show()