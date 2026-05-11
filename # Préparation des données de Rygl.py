# Préparation des données de Rygl

from astropy.coordinates import SkyCoord, Galactic
import astropy.units as u

ra_rygl = 307.711 # degrés
dec_rygl = 41.041 # degrés
mu_alpha_rygl = -2.84 # mas/yr (attention, c'est souvent mu_alpha * cos(delta) dans ces articles)
mu_delta_rygl = -4.14 # mas/yr

# 1. On crée la coordonnée dans le repère Équatorial (ICRS)
coord_eq = SkyCoord(
    ra=ra_rygl * u.deg, 
    dec=dec_rygl * u.deg,
    pm_ra_cosdec=mu_alpha_rygl * (u.mas / u.yr),
    pm_dec=mu_delta_rygl * (u.mas / u.yr),      
    frame='icrs'
)

# On transforme le tout en repère Galactique
coord_gal = coord_eq.transform_to(Galactic)

# On extrait les vraies valeurs mu_l et mu_b
mu_l_converti = coord_gal.pm_l_cosb.value
mu_b_converti = coord_gal.pm_b.value

print(f"mu_l = {mu_l_converti:.2f} mas/yr")
print(f"mu_b = {mu_b_converti:.2f} mas/yr")