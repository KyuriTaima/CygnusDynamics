# Graph creation of proper motion

import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
import pandas as pd

# Import csv data
objects_data = pd.read_csv('Cygnus_Objects_Datas_Uncertainties_v3.csv')

# Extract relevant columns
names = objects_data['Object_Name'].values
mu_l = objects_data['mu_l_mas_yr'].values
mu_b = objects_data['mu_b_mas_yr'].values
l_deg = objects_data['l_deg'].values

# Color assignation
point_colors = ['red' if 'Group' in str(nom) else 'blue' for nom in names]

# Plot mu_b as a function of mu_l
plt.figure(figsize=(10, 6))
plt.scatter(mu_l, mu_b, c=point_colors)
plt.scatter([], [], color='red', label="OB Stars Clusters")
plt.scatter([], [], color='blue', label="Molecular Clouds")
plt.xlabel('Proper Motion in Galactic Longitude ($\mu_l$) [mas/yr]')
plt.ylabel('Proper Motion in Galactic Latitude ($\mu_b$) [mas/yr]')
plt.title('Proper Motion in Galactic Coordinates')

# Add labels for each point with an offset
for i in range(len(names)):
    plt.annotate(
        names[i], 
        (mu_l[i], mu_b[i]), 
        textcoords="offset points", 
        xytext=(5, 5),              
        fontsize=8, 
        ha='left',                   
        color=point_colors[i]        
    )

plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()
plt.show()


# Plot mu_l as a function of l
plt.figure(figsize=(10, 6))
plt.scatter(l_deg, mu_l, c=point_colors)
plt.scatter([], [], color='red', label="OB Stars Clusters")
plt.scatter([], [], color='blue', label="Molecular Clouds")
plt.xlabel('Galactic Longitude ($l$) [degrees]')
plt.ylabel('Proper Motion in Galactic Longitude ($\mu_l$) [mas/yr]')
plt.title('Proper Motion in Galactic Longitude vs Galactic Longitude')

# Add labels for each point with an offset
for i in range(len(names)):
    plt.annotate(
        names[i], 
        (l_deg[i], mu_l[i]), 
        textcoords="offset points", 
        xytext=(5, 5),              
        fontsize=8, 
        ha='left', 
        color=point_colors[i]
    )
plt.grid(True, linestyle='--', alpha=0.6)
plt.xlim(83, 71)
plt.legend()
plt.show()