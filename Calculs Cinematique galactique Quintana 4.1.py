# Changement de référentiel des vitesses des groupes OB de la région du Cygne

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord, Galactic

# Intialisation des paramètres de la simulation
R0 = 8.15        # kpc Distance du soleil au centre galactique
theta_helio = 229.0 # km/s Vitesse de rotation des objets au niveau du soleil autour du centre galactique (Eilers et Al. 2019)
l_deg = 72 # degrés de longitude galactique de l'amas étudié
d = 1.8 # kpc Distance de l'amas étudié au soleil
angle_visee_intercept_deg = 4.5 # Angle entre la vitesse de l'objet et la ligne de visée, en degrés

# Schönrich et al. (2010) : Vitesse du Soleil par rapport au LSR
U_sun = 11.10
V_sun = 12.24
W_sun = 7

# Calculs préliminaires
l_rad = np.radians(l_deg) # Conversion de la longitude galactique en radians
angle_visee_intercept_rad = np.radians(angle_visee_intercept_deg) # Conversion de l'angle de visée en radians

def rotation_curve(R, R0, Theta0):
    #Modèle de courbe de rotation galactique utilisé par Quintana et Al. (2021)
    v = Theta0 - 1.7*(R-R0)
    return v

# Calcul de la vitesse de rotation à la position de l'amas
R = np.sqrt(R0**2 + d**2 - 2 * R0 * d * np.cos(l_rad)) # Distance de l'amas au centre galactique en utilisant le théorème d'Al-Kashi
print(f"Distance de l'amas au centre galactique : {R:.2f} kpc")
Theta_R = rotation_curve(R, R0, theta_helio) # Vitesse de rotation à la position de l'amas
print(f"Vitesse de rotation à la position de l'amas : {Theta_R:.2f} km/s")

v_trans_galactique = - (Theta_R * np.sin(angle_visee_intercept_rad))
print(f"Vitesse transversale (rotation) : {v_trans_galactique:.2f} km/s")

# Correction de la vitesse du soleil sur le plan de vision de l'amas
v_sun_projection = U_sun * np.sin(l_rad) - V_sun * np.cos(l_rad)
print(f"Vitesse de projection du soleil : {v_sun_projection:.2f} km/s")

# On retire la vitesse du soleil pour obtenir la vitesse héliocentrique brute de l'amas
v_helio = v_trans_galactique - v_sun_projection
print(f"Vitesse héliocentrique brute de l'amas : {v_helio:.2f} km/s")

# Transformation en mouvement propre
mu_l = v_helio / (4.74 * d) # Conversion de la vitesse héliocentrique en mouvement propre en mas/yr
print(f"Mouvement propre en longitude galactique : {mu_l:.2f} mas/yr")

# Mouvement propre en latitude
V_trans_b = - W_sun
mu_b = V_trans_b / (4.74 * d) # Conversion de la vitesse transversale en latitude en mouvement propre en mas/yr
print(f"Mouvement propre en latitude galactique : {mu_b:.2f} mas/yr")


