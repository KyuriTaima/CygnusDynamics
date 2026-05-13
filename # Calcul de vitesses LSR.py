import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord, Galactic

# Initialisation des paramètres
R0 = 8.277        # kpc Distance du soleil au centre galactique
Theta0 = 234.8  # km/s Vitesse de rotation fixe des objets autour du centre galactique

# Schönrich et al. (2010) : Vitesse du Soleil par rapport au LSR
U_sun = 9.1 * u.km / u.s
V_sun = 8.66 * u.km / u.s
W_sun = 7.90 * u.km / u.s

# tableau de valeurs pour les différentes sources dans Cygnus X: W75N, DR21, DR20, IRAS 20290+4052

ra = [309.651, 309.758, 309.254, 307.711] * u.deg
dec = [42.626, 42.416, 41.582, 41.041] * u.deg
d = [1.30, 1.50, 1.46, 1.36]
pm_ra = [-1.97, -2.84, -3.29, -2.84] * u.mas / u.yr
pm_dec = [-4.16, -3.80, -4.83, -4.14] * u.mas / u.yr
v_lsr_radio = [9.0, -3.0, -3.0, -1.4] * u.km / u.s
l = [81.87, 81.75, 80.86, 79.74]
b = [0.78, 0.59, 0.38, 0.99]

# Note sur le v_lsr : 
# Les données radio v_lsr utilisent souvent le "Standard Solar Motion" historique 
# (U=10.3, V=15.3, W=7.7) pour définir le LSR d'origine dans leurs télescopes.
# Pour retrouver le V_helio brut, il faut utiliser cette ancienne norme.
U_old = 10.3 * u.km / u.s
V_old = 15.3 * u.km / u.s
W_old = 7.7 * u.km / u.s

for i in range(len(ra)):
    print(f"\n--- Source {i+1} ---")
    print(f"RA: {ra[i]}, Dec: {dec[i]}, Distance: {d[i]}, PM_RA: {pm_ra[i]}, PM_Dec: {pm_dec[i]}, V_LSR_radio: {v_lsr_radio[i]}")
    coord_base = SkyCoord(ra=ra[i], dec=dec[i], frame='icrs')
    gal_base = coord_base.transform_to(Galactic)
    l_rad = gal_base.l.radian
    b_rad = gal_base.b.radian

    # Récupération de la vitesse héliocentrique brute (avec l'ancien Soleil radio)
    # c'est la vitesse de l'objet, par rapport au point LSR, en ayant corrigé de la vitesse du soleil
    v_sun_proj_old = (U_old * np.cos(l_rad) * np.cos(b_rad) + 
                    V_old * np.sin(l_rad) * np.cos(b_rad) + 
                    W_old * np.sin(b_rad))
    # En soustrayant cette projection de la vitesse LSR radio, on obtient la vitesse héliocentrique brute de l'objet, sans les corrections du NOUVEAU Soleil.
    v_helio = v_lsr_radio[i] - v_sun_proj_old

    # Utilisation d'Astropy pour les composantes Héliocentriques
    c_complet = SkyCoord(ra=ra[i], dec=dec[i], distance=d[i]*u.kpc,
                        pm_ra_cosdec=pm_ra[i], pm_dec=pm_dec[i],
                        radial_velocity=v_helio, frame='icrs')
    gal = c_complet.transform_to(Galactic)

    # Application des composantes du Soleil (Schönrich 2010)
    U_lsr_sun = gal.velocity.d_x.value + U_sun.value
    V_lsr_sun = gal.velocity.d_y.value + V_sun.value
    W_lsr_sun = gal.velocity.d_z.value + W_sun.value

    # La transformation Galactocentrique (Reid et al. methodology)
    # Passage au repère galactocentrique fixe
    U_gal = U_lsr_sun # meme direction (radiale vers le centre galactique), donc pas de changement
    V_gal = V_lsr_sun + Theta0 # On ajoute la vitesse de rotation du LSR pour passer à un repère galactocentrique fixe

    # Calcul de l'angle Beta
    X = R0 - d[i] * np.cos(l_rad) * np.cos(b_rad)
    Y = d[i] * np.sin(l_rad) * np.cos(b_rad)
    beta = np.arctan2(Y, X) # Angle de rotation

    # Rotation des axes et soustraction de la rotation homogène (Theta0)
    U_pec = U_gal * np.cos(beta) - V_gal * np.sin(beta)
    V_pec = U_gal * np.sin(beta) + V_gal * np.cos(beta) - Theta0
    W_pec = W_lsr_sun

    # Calcul des composantes de vitesses longitudinales et latitudales des objets
    v_l = -U_pec*np.sin(np.radians(l[i]))+V_pec*np.cos(np.radians(l[i]))
    v_b = -U_pec*np.cos(np.radians(l[i]))*np.sin(np.radians(b[i]))-V_pec*np.sin(np.radians(l[i]))*np.sin(np.radians(b[i]))+W_pec*np.cos(np.radians(b[i]))
    v_norm = np.sqrt(v_l**2 + v_b**2)
    # Calcul de leur norme en cm pour placement sur la figure 12 de Rygl et Al. 2018
    v_l_cm = (2.3/10) * v_l
    v_b_cm = (2.3/10) * v_b
    v_cm = (2.3/10) * v_norm

    print(f"--- RÉSULTATS AVEC R0={R0}, THETA0={Theta0} ---")
    print(f"U = {U_pec:.2f} km/s")
    print(f"V = {V_pec:.2f} km/s")
    print(f"W = {W_pec:.2f} km/s")
    print(f"v_helio = {v_helio:.2f} km/s")
    print(f"Vitesses sur la figure 12: v_l = {v_l:.2f} km/s; v_b = {v_b:.2f} km/s; v_cm = {v_cm:.2f} cm")