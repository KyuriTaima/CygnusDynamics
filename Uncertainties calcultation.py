# Uncertainties calcultation

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord, Galactic
import pandas as pd

# Initialisation des paramètres

r0 = 8.15        # kpc Distance du soleil au centre galactique
theta0 = 236.0  # km/s Vitesse de rotation fixe des objets autour du centre galactique
theta0_err = 7.0 # km/s Incertitude sur la vitesse de rotation homogène de la galaxie au LSR
# Schönrich et al. (2010) : Vitesse du Soleil par rapport au LSR
U_sun = 11.1 # km/s
V_sun = 15.0 # km/s
W_sun = 7.2 # km/s
U_sun_err = 1.2 # km/s
V_sun_err = 10.0 # km/s
W_sun_err = 1.1 # km/s

# Ancienne norme de vitesse du soleil utilisée pour les données radio v_lsr
U_old = 10.3 # km/s
V_old = 15.3 # km/s
W_old = 7.7 # km/s

names = ["Group A", "Group B", "Group C", "Group D", "Group E", "Group F", "W75N", "DR21", "DR20", "IRAS20290+4052"]

# --- Coordonnées et distances ---
right_ascension = [301.45, 304.37, 305.47, 304.34, 308.08, 302.95, 309.651, 309.758, 309.254, 307.711] # degrés
declination = [35.68, 41.43, 37.87, 37.64, 41.30, 36.58, 42.626, 42.416, 41.582, 41.041] # degrés
galactic_longitude = [72.61, 78.58, 76.11, 75.44, 80.19, 74.04, 81.87, 81.75, 80.86, 79.74] # degrés
galactic_latitude = [2.06, 3.31, 0.54, 1.19, 0.85, 1.44, 0.78, 0.59, 0.38, 0.99] # degrés
distance = [1894.5, 1726.3, 1713.1, 2000.1, 1674.0, 1985.2, 1300, 1500, 1460, 1360] # pc
distance_err = [4, 7, 10, 10, 10, 10 , 0.07, 0.08, 0.09, 0.12] # pc

# --- Mouvements propres ---
mu_l = [-6.90, -5.47, -6.57, -5.55, -4.72, -6.11, -4.50, -4.74, -5.84, -5.02] # mas/yr
mu_l_err = [0.24, 0.34, 0.27, 0.16, 0.27, 0.30, 0.01, 0.01, 0.01, 0.02] # mas/yr

mu_b = [-1.35, -0.59, -1.19, -1.34, -0.96, -1.33, -0.96, -0.06, -0.29, -0.15] # mas/yr
mu_b_err = [0.17, 0.27, 0.12, 0.14, 0.25, 0.10, 0.01, 0.01, 0.01, 0.02] # mas/yr

# --- Dispersions de vitesse ---
sigma_v_l = [2.20, 2.72, 2.24, 1.53, 2.15, 2.80, np.nan, np.nan, np.nan, np.nan] # km/s
sigma_v_l_err_up = [0.55, 0.47, 0.46, 0.45, 0.66, 0.58, np.nan, np.nan, np.nan, np.nan] # km/s
sigma_v_l_err_down = [0.47, 0.37, 0.35, 0.41, 0.45, 0.48, np.nan, np.nan, np.nan, np.nan]  # km/s

sigma_v_b = [1.55, 2.18, 0.96, 1.32, 1.96, 1.50, np.nan, np.nan, np.nan, np.nan]  # km/s
sigma_v_b_err_up = [0.32, 0.41, 0.22, 0.30, 0.52, 0.30, np.nan, np.nan, np.nan, np.nan] # km/s
sigma_v_b_err_down = [0.24, 0.33, 0.18, 0.26, 0.29, 0.25, np.nan, np.nan, np.nan, np.nan] # km/s

# --- Comptages d'étoiles ---
observed_b_stars = [112, 82, 76, 75, 110, 139, np.nan, np.nan, np.nan, np.nan]
observed_b_stars_err_up = [3, 3, 2, 2, 4, 3, np.nan, np.nan, np.nan, np.nan]
observed_b_stars_err_down = [3, 2, 3, 3, 3, 3, np.nan, np.nan, np.nan, np.nan]

observed_o_stars = [3, 4, 10, 4, 35, 5, np.nan, np.nan, np.nan, np.nan]
observed_o_stars_err_up = [1, 1, 1, 2, 2, 2, np.nan, np.nan, np.nan, np.nan]
observed_o_stars_err_down = [1, 1, 2, 1, 3, 1, np.nan, np.nan, np.nan, np.nan]

corrected_o_stars = [13, 11, 12, 11, 23, 17, np.nan, np.nan, np.nan, np.nan]
corrected_o_stars_err_up = [4, 4, 4, 4, 6, 4, np.nan, np.nan, np.nan, np.nan]
corrected_o_stars_err_down = [3, 3, 3, 3, 5, 5, np.nan, np.nan, np.nan, np.nan]

# --- Masse totale estimée ---
total_stellar_mass = [2344, 1978, 2165, 1569, 4198, 2955, np.nan, np.nan, np.nan, np.nan] # Msun
total_stellar_mass_err_up = [275, 245, 245, 215, 354, 349, np.nan, np.nan, np.nan, np.nan] # Msun
total_stellar_mass_err_down = [252, 224, 224, 192, 319, 269, np.nan, np.nan, np.nan, np.nan] # Msun

# --- Gradients de vitesse ---
# Gradient selon (l) en mas/yr/deg
vel_grad_l_mas = [0.24, 0.39, 0.38, 0.20, 0.43, 0.35, np.nan, np.nan, np.nan, np.nan]
vel_grad_l_mas_err = [0.07, 0.04, 0.04, 0.05, 0.06, 0.03, np.nan, np.nan, np.nan, np.nan]

# Gradient selon (l) en km/s/pc
vel_grad_l_km = [0.07, 0.11, 0.10, 0.05, 0.12, 0.09, np.nan, np.nan, np.nan, np.nan]
vel_grad_l_km_err = [0.02, 0.01, 0.01, 0.01, 0.01, 0.01, np.nan, np.nan, np.nan, np.nan]

# Gradient selon (b) en mas/yr/deg
vel_grad_b_mas = [0.34, -0.03, 0.10, 0.07, 0.16, 0.06, np.nan, np.nan, np.nan, np.nan]
vel_grad_b_mas_err = [0.04, 0.04, 0.04, 0.03, 0.05, 0.03, np.nan, np.nan, np.nan, np.nan]

# Gradient selon (b) en km/s/pc
vel_grad_b_km = [0.09, -0.01, 0.03, 0.02, 0.04, 0.02, np.nan, np.nan, np.nan, np.nan]
vel_grad_b_km_err = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, np.nan, np.nan, np.nan, np.nan]

# --- Âge d'expansion ---
expansion_age_l = [13.98, 8.89, 7.93, 19.57, 8.37, 10.87, np.nan, np.nan, np.nan, np.nan] # Myr
expansion_age_l_err = [3.99, 0.81, 0.80, 3.91, 0.70, 1.21, np.nan, np.nan, np.nan, np.nan] # Myr

expansion_age_b = [10.87, np.nan, 32.62, np.nan, 24.46, np.nan, np.nan, np.nan, np.nan, np.nan] # Myr
expansion_age_b_err = [1.21, np.nan, 10.87, np.nan, 6.11, np.nan, np.nan, np.nan, np.nan, np.nan] # Myr

# --- Vitesses Radiales (Ligne de visée / Héliocentriques) ---
# Données issues de la Table 2 de Quintana & Wright (2022)

radial_velocity = [-12.50, -21.00, -19.20, -7.00, -12.35, -13.48, -7.71, -19.70, -19.79, -18.42] # km/s
radial_velocity_err = [2.77, 10.80, 2.98, 3.88, 2.06, 5.06, 0.01, 0.01, 0.01, 0.01] # km/s

v_lsr_err_list = [2.77, 10.80, 2.98, 3.88, 2.06, 5.06, 0.01, 0.01, 0.01, 0.02] # km/s

print(f"--- PARAMÈTRES INITIAUX ---")
print(f"R0 = {r0} kpc")
print(f"Theta0 = {theta0} km/s")
print(f"U_sun = {U_sun} km/s")
print(f"V_sun = {V_sun} km/s")
print(f"W_sun = {W_sun} km/s")
print(f"----------------------------")
print(f"----------------------------")
print(f"----------------------------")

# Initialisation des listes pour stocker les résultats des vitesses dans les différents référentiels afin de pouvoir les exporter ensuite
U_obs_list = []
V_obs_list = []
W_obs_list = []
U_lsr_list = []
V_lsr_list = []
W_lsr_list = []
U_pec_list = []
V_pec_list = []
W_pec_list = []
U_lsr_abs_list = [] # Vitesse lsr absolue en rajoutant la vitesse de rotation homogène theta0
V_lsr_abs_list = []
W_lsr_abs_list = []
beta_list = []
V_gal_abs_list = []
v_lsr_list = [] # Vitesses lsr comparables avec les données radio v_lsr de Rygl et Al. 2018 tableau 2
v_l_pec_kms_list = [] # Vitesse longitudinale
v_b_pec_kms_list = [] # Vitesse latitudinale
v_lb_kms_list = [] #norme de la vitesse dans le plan galactique
# Initialisation des coordonnées spatiales héliocentriques
X_helio_list = []
Y_helio_list = []
Z_helio_list = []
X_helio_err_list = []
Y_helio_err_list = []
Z_helio_err_list = []
U_obs_err_list = []
V_obs_err_list = []
W_obs_err_list = []
U_lsr_err_list = []
V_lsr_err_list = []
W_lsr_err_list = []
U_lsr_abs_err_list = []
V_lsr_abs_err_list = []
W_lsr_abs_err_list = []
U_pec_err_list = []
V_pec_err_list = []
W_pec_err_list = []
V_gal_abs_err_list = [] # Liste pour stocker les incertitudes sur les vitesses absolues dans le référentiel galactocentrique
v_l_pec_err_list = [] # Liste pour stocker les incertitudes sur les vitesses longitudinales
v_b_pec_err_list = [] # Liste pour stocker les incertitudes sur les vitesses latitudinales
v_lb_err_list = [] # Liste pour stocker les incertitudes sur les normes de vitesses dans le plan galactique

# Boucle sur les groupes pour calculer les composantes de vitesses dans le référentiel LSR et Galactocentrique
for i in range(len(names)):
    print(f"\n--- {names[i]} ---")
    # Création de l'objet Astropy dans le référentiel Galactique
    coord = SkyCoord(
            frame=Galactic,
            l=galactic_longitude[i] * u.deg, b=galactic_latitude[i] * u.deg,
            distance=distance[i] * u.pc,
            pm_l_cosb=mu_l[i] * (u.mas / u.yr), 
            pm_b=mu_b[i] * (u.mas / u.yr),
            radial_velocity=radial_velocity[i] * (u.km / u.s)
        )
    
    U_obs = coord.velocity.d_x.value 
    V_obs = coord.velocity.d_y.value
    W_obs = coord.velocity.d_z.value
    U_obs_list.append(U_obs)
    V_obs_list.append(V_obs)
    W_obs_list.append(W_obs)

    X_helio = coord.cartesian.x.value
    Y_helio = coord.cartesian.y.value
    Z_helio = coord.cartesian.z.value
    
    X_helio_list.append(X_helio)
    Y_helio_list.append(Y_helio)
    Z_helio_list.append(Z_helio)

    # Correction de la vitesse du soleil pour obtenir les composantes dans le référentiel LSR
    U_lsr = U_obs + U_sun
    V_lsr = V_obs + V_sun
    W_lsr = W_obs + W_sun
    U_lsr_list.append(U_lsr)
    V_lsr_list.append(V_lsr)
    W_lsr_list.append(W_lsr)

    # Calcul de l'angle Beta pour la transformation en référentiel galactocentrique
    l_rad = np.radians(galactic_longitude[i])
    b_rad = np.radians(galactic_latitude[i])
    d_val = distance[i] * (u.pc.to(u.kpc)) # conversion de la distance en kpc pour le calcul de l'angle Beta

    X = r0 - d_val * np.cos(l_rad) * np.cos(b_rad)
    Y = d_val * np.sin(l_rad) * np.cos(b_rad)
    beta = np.arctan2(Y, X) # Angle de rotation
    beta_list.append(beta)

    # Vitesse lsr absolue en rajoutant la vitesse de rotation homogène theta0
    U_lsr_abs = U_lsr # pas changement pour la composante radiale U
    V_lsr_abs = V_lsr + theta0 # On ajoute la vitesse de rotation du LSR pour passer à un repère galactocentrique fixe
    W_lsr_abs = W_lsr # pas de changement pour la composante verticale W
    U_lsr_abs_list.append(U_lsr_abs)
    V_lsr_abs_list.append(V_lsr_abs)
    W_lsr_abs_list.append(W_lsr_abs)

    # Rotation des axes pour obtenir les composantes de vitesses dans le référentiel galactocentrique
    U_pec = U_lsr_abs * np.cos(beta) - V_lsr_abs * np.sin(beta)
    V_pec = U_lsr_abs * np.sin(beta) + V_lsr_abs * np.cos(beta) - theta0
    W_pec = W_lsr_abs # pas de changement pour la composante verticale W
    U_pec_list.append(U_pec)
    V_pec_list.append(V_pec)
    W_pec_list.append(W_pec)

    V_gal_abs = V_pec + theta0 # Vitesse de rotation homogène à rajouter pour obtenir la vitesse absolue dans le référentiel galactocentrique
    V_gal_abs_list.append(V_gal_abs)

    print(f"Vitesses dans le référentiel LSR : U_lsr = {U_lsr:.2f} km/s, V_lsr = {V_lsr:.2f} km/s, W_lsr = {W_lsr:.2f} km/s")
    print(f"Vitesses dans le référentiel galactocentrique : U_pec = {U_pec:.2f} km/s, V_pec = {V_pec:.2f} km/s, W_pec = {W_pec:.2f} km/s")

    # Calcul de la vitesse v_lsr comparable avec l'article de Rygl et Al. 2018 tableau 2
    v_sun_proj_old = (U_old * np.cos(l_rad) * np.cos(b_rad) + 
                      V_old * np.sin(l_rad) * np.cos(b_rad) + 
                      W_old * np.sin(b_rad))
    v_lsr = radial_velocity[i] + v_sun_proj_old
    v_lsr_list.append(v_lsr)

    # Calcul des vitesses longitudinales et latitudinales
    v_l_pec = -U_pec * np.sin(l_rad) + V_pec * np.cos(l_rad)
    v_b_pec = (-U_pec * np.cos(l_rad) * np.sin(b_rad) 
               - V_pec * np.sin(l_rad) * np.sin(b_rad) 
               + W_pec * np.cos(b_rad))
    
    v_l_pec_kms_list.append(v_l_pec)
    v_b_pec_kms_list.append(v_b_pec)
    v_lb_kms = np.sqrt(v_l_pec**2 + v_b_pec**2)
    v_lb_kms_list.append(v_lb_kms)

# Calcul des incertitudes
# Object definition
n_samples = 10000

for i in range(len(names)):
    print(f"\n--- {names[i]} ---")
    distances = np.random.normal(distance[i], distance_err[i], n_samples) * u.pc
    longitudes = np.random.normal(galactic_longitude[i], 0.01, n_samples) * u.deg
    latitudes = np.random.normal(galactic_latitude[i], 0.01, n_samples) * u.deg

    # Feed the arrays to SkyCoord
    c = SkyCoord(frame=Galactic, l=longitudes, b=latitudes, distance=distances)

    # Transform everything at once
    cart = c.cartesian

    # Result
    X_helio_err = cart.x.std()
    Y_helio_err = cart.y.std()
    Z_helio_err = cart.z.std()

    X_helio_err_list.append(X_helio_err)
    Y_helio_err_list.append(Y_helio_err)
    Z_helio_err_list.append(Z_helio_err)

    # Monte-Carlo simulation for heliocentric velocity components uncertainties
    mu_l_samples = np.random.normal(mu_l[i], mu_l_err[i], n_samples) * (u.mas / u.yr)
    mu_b_samples = np.random.normal(mu_b[i], mu_b_err[i], n_samples) * (u.mas / u.yr)
    radial_velocity_samples = np.random.normal(radial_velocity[i], radial_velocity_err[i], n_samples) * (u.km / u.s)
    c_samples = SkyCoord(frame=Galactic, l=longitudes, b=latitudes, distance=distances,
                         pm_l_cosb=mu_l_samples, pm_b=mu_b_samples, radial_velocity=radial_velocity_samples)
    U_obs_samples = c_samples.velocity.d_x.value
    V_obs_samples = c_samples.velocity.d_y.value
    W_obs_samples = c_samples.velocity.d_z.value
    U_obs_err = np.std(U_obs_samples)
    V_obs_err = np.std(V_obs_samples)
    W_obs_err = np.std(W_obs_samples)
    U_obs_err_list.append(U_obs_err)
    V_obs_err_list.append(V_obs_err)
    W_obs_err_list.append(W_obs_err)

    #Taylor expansion for LSR velocity components uncertainties
    U_lsr_err = np.sqrt(U_obs_err**2 + U_sun_err**2)
    V_lsr_err = np.sqrt(V_obs_err**2 + V_sun_err**2)
    W_lsr_err = np.sqrt(W_obs_err**2 + W_sun_err**2)
    U_lsr_err_list.append(U_lsr_err)
    V_lsr_err_list.append(V_lsr_err)
    W_lsr_err_list.append(W_lsr_err)

    # Taylor expansion for LSR non corrected by the galactic velocity
    U_lsr_abs_err = U_lsr_err
    V_lsr_abs_err = np.sqrt(V_lsr_err**2 + theta0_err**2)
    W_lsr_abs_err = W_lsr_err
    U_lsr_abs_err_list.append(U_lsr_abs_err)
    V_lsr_abs_err_list.append(V_lsr_abs_err)
    W_lsr_abs_err_list.append(W_lsr_abs_err)

    # Taylor expansion for peculiar velocity components uncertainties
    beta = beta_list[i]
    U_pec_err = np.sqrt((U_lsr_abs_err * np.cos(beta))**2 + (V_lsr_abs_err * np.sin(beta))**2)
    V_pec_err = np.sqrt((U_lsr_abs_err * np.sin(beta))**2 + (V_lsr_abs_err * np.cos(beta))**2 + theta0_err**2)
    W_pec_err = W_lsr_abs_err
    U_pec_err_list.append(U_pec_err)
    V_pec_err_list.append(V_pec_err)
    W_pec_err_list.append(W_pec_err)

    # Taylor expansion for galactocentric absolute velocity
    V_gal_abs_err = np.sqrt(V_pec_err**2 + theta0_err**2)
    V_gal_abs_err_list.append(V_gal_abs_err)

    # Taylor expansion for longitudinal and latitudinal peculiar velocities
    v_l_pec_err = np.sqrt((U_pec_err * np.sin(l_rad))**2 + (V_pec_err * np.cos(l_rad))**2)
    v_b_pec_err = np.sqrt((U_pec_err * np.cos(l_rad) * np.sin(b_rad))**2 + (V_pec_err * np.sin(l_rad) * np.sin(b_rad))**2 + (W_pec_err * np.cos(b_rad))**2)
    v_l_pec_err_list.append(v_l_pec_err)
    v_b_pec_err_list.append(v_b_pec_err)
    v_lb_err = np.sqrt(v_l_pec_kms_list[i]**2 * v_l_pec_err**2 + v_b_pec_kms_list[i]**2 * v_b_pec_err**2) / v_lb_kms_list[i]
    v_lb_err_list.append(v_lb_err)
    

# --- CRÉATION DU DATAFRAME PANDAS ---
# Mapping de chaque colonne avec les données correspondantes dans le script
data = {
    "Object_Name": names,
    "RA_deg": right_ascension,
    "Dec_deg": declination,
    "l_deg": galactic_longitude,
    "b_deg": galactic_latitude,
    "X_helio_pc": X_helio_list,
    "Y_helio_pc": Y_helio_list,
    "Z_helio_pc": Z_helio_list,
    "Distance_pc": distance,
    "mu_l_mas_yr": mu_l,
    "mu_b_mas_yr": mu_b,
    "Radial_Velocity_kms": radial_velocity,
    "Radial_Velocity_err": radial_velocity_err,
    "Sigma_v_l_kms": sigma_v_l,
    "Sigma_v_l_err_up": sigma_v_l_err_up,
    "Sigma_v_l_err_down": sigma_v_l_err_down,
    "Sigma_v_b_kms": sigma_v_b,
    "Sigma_v_b_err_up": sigma_v_b_err_up,
    "Sigma_v_b_err_down": sigma_v_b_err_down,
    "N_B_stars": observed_b_stars,
    "N_B_stars_err_up": observed_b_stars_err_up,
    "N_B_stars_err_down": observed_b_stars_err_down,
    "N_O_stars": observed_o_stars,
    "N_O_stars_err_up": observed_o_stars_err_up,
    "N_O_stars_err_down": observed_o_stars_err_down,
    "N_O_stars_corrected": corrected_o_stars,
    "N_O_stars_corrected_err_up": corrected_o_stars_err_up,
    "N_O_stars_corrected_err_down": corrected_o_stars_err_down,
    "Total_Stellar_Mass_Msun": total_stellar_mass,
    "Total_Stellar_Mass_err_up": total_stellar_mass_err_up,
    "Total_Stellar_Mass_err_down": total_stellar_mass_err_down,
    "Vel_Grad_l_mas_deg": vel_grad_l_mas,
    "Vel_Grad_l_mas_deg_err": vel_grad_l_mas_err,
    "Vel_Grad_l_km_pc": vel_grad_l_km,
    "Vel_Grad_l_km_pc_err": vel_grad_l_km_err,
    "Vel_Grad_b_mas_deg": vel_grad_b_mas,
    "Vel_Grad_b_mas_deg_err": vel_grad_b_mas_err,
    "Vel_Grad_b_km_pc": vel_grad_b_km,
    "Vel_Grad_b_km_pc_err": vel_grad_b_km_err,
    "Expansion_Age_l_Myr": expansion_age_l,
    "Expansion_Age_l_Myr_err": expansion_age_l_err,
    "Expansion_Age_b_Myr": expansion_age_b,
    "Expansion_Age_b_Myr_err": expansion_age_b_err,
    "U_obs_kms": U_obs_list,
    "V_obs_kms": V_obs_list,
    "W_obs_kms": W_obs_list,
    "U_lsr_kms": U_lsr_list,
    "V_lsr_kms": V_lsr_list,
    "W_lsr_kms": W_lsr_list,
    "U_pec_kms": U_pec_list,
    "V_pec_kms": V_pec_list,
    "W_pec_kms": W_pec_list,
    "U_lsr_abs_kms": U_lsr_abs_list,
    "V_lsr_abs_kms": V_lsr_abs_list,
    "W_lsr_abs_kms": W_lsr_abs_list,
    "Beta_rad": beta_list,
    "v_lsr_kms": v_lsr_list,
    "V_gal_abs_kms": V_gal_abs_list,
    "v_l_kms": v_l_pec_kms_list,
    "v_b_kms": v_b_pec_kms_list,
    "v_lb_kms": v_lb_kms_list,
    "Distance_err_pc": distance_err,
    "RA_err_deg": [0.00001] * len(names),
    "Dec_err_deg": [0.00001] * len(names),
    "l_err_deg": [0.00001] * len(names),
    "b_err_deg": [0.00001] * len(names), 
    "X_helio_err_pc": X_helio_err_list,
    "Y_helio_err_pc": Y_helio_err_list,
    "Z_helio_err_pc": Z_helio_err_list,
    "mu_l_err": mu_l_err,
    "mu_b_err": mu_b_err,
    "Radial_Velocity_err_kms": radial_velocity_err,
    "U_obs_err_kms": U_obs_err_list,
    "V_obs_err_kms": V_obs_err_list,
    "W_obs_err_kms": W_obs_err_list,
    "U_lsr_err_kms": U_lsr_err_list,
    "V_lsr_err_kms": V_lsr_err_list,
    "W_lsr_err_kms": W_lsr_err_list,
    "U_lsr_abs_err_kms": U_lsr_abs_err_list,
    "V_lsr_abs_err_kms": V_lsr_abs_err_list,
    "W_lsr_abs_err_kms": W_lsr_abs_err_list,
    "U_pec_err_kms": U_pec_err_list,
    "V_pec_err_kms": V_pec_err_list,
    "W_pec_err_kms": W_pec_err_list,
    "v_lsr_err_kms": v_lsr_err_list,
    "V_gal_abs_err_kms": V_gal_abs_err_list,
    "v_l_pec_err_kms": v_l_pec_err_list,
    "v_b_pec_err_kms": v_b_pec_err_list,
    "v_lb_err_kms": v_lb_err_list
}

# Display the length of each list to ensure they match
for key, value in data.items():
    print(f"{key}: {len(value)}")

df = pd.DataFrame(data)

# Export en CSV
csv_filename = "Cygnus_Objects_Datas_Uncertainties.csv"
df.to_csv(csv_filename, index=False, float_format='%.2f')

print(f"\nDataFrame créé et exporté en CSV : {csv_filename}")

