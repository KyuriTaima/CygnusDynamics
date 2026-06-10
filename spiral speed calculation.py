import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress

# --- Data Loading and Initialization ---
df = pd.read_csv('Cygnus_Objects_Datas_Uncertainties.csv')
df = df.dropna(subset=['Distance_pc', 'vx_gc_kms', 'vy_gc_kms'])

object_names = df['Object_Name'].values
distances_pc = df['Distance_pc'].values
l_rad = np.radians(df['l_deg'].values)
b_rad = np.radians(df['b_deg'].values)

# Absolute Galactocentric velocities from CSV
vx_gal = df['vx_gc_kms'].values
vy_gal = df['vy_gc_kms'].values

# Calculate current global Galactocentric coordinates (pc)
R0 = 8277 
x_gal = R0 - distances_pc * np.cos(l_rad) * np.cos(b_rad)
y_gal = distances_pc * np.sin(l_rad) * np.cos(b_rad)

# --- Define Chronological Generations ---
generations = {
    0.0: ["W75N", "DR21", "DR20", "IRAS20290+4052"],
    4.5: ["Group B", "Group C", "Group E"],
    8.5: ["Group A", "Group D", "Group F"]
}

conversion_pc_myr = 1.022
lookback_times = []
mean_azimuths = []

# --- Global Traceback Loop ---
for t_lookback, names in generations.items():
    indices = [list(object_names).index(name) for name in names if name in object_names]
    
    # Trace back positions in the absolute Galactocentric frame
    x_past = x_gal[indices] - vx_gal[indices] * t_lookback * conversion_pc_myr
    y_past = y_gal[indices] - vy_gal[indices] * t_lookback * conversion_pc_myr
    
    # Calculate absolute Galactic azimuth (rad)
    phi_past = np.arctan2(y_past, x_past)
    
    lookback_times.append(t_lookback)
    mean_azimuths.append(np.mean(phi_past))

# --- Spiral Pattern Speed Calculation ---
# Linear regression: phi = intercept + slope * t_lookback
slope, intercept, r_value, p_value, std_err = linregress(lookback_times, mean_azimuths)

# Omega_p = -slope. Convert rad/Myr to km/s/kpc
rad_myr_to_kms_kpc = 977.81
omega_p = -slope * rad_myr_to_kms_kpc
omega_p_err = std_err * rad_myr_to_kms_kpc

print(f"--- Spiral Wave Propagation Results ---")
print(f"Spiral Pattern Speed (Omega_p): {omega_p:.2f} +/- {omega_p_err:.2f} km/s/kpc")
print(f"R-squared value               : {r_value**2:.4f}")
print("-" * 40)

# --- Verification Plot ---
plt.figure(figsize=(8, 5))
plt.scatter(lookback_times, mean_azimuths, color='crimson', s=100, zorder=3, label='Observed Generations')
plt.plot(lookback_times, intercept + slope * np.array(lookback_times), color='royalblue', lw=2, label='Linear Fit')

plt.title('Galactocentric Azimuthal Drift Over Lookback Time')
plt.xlabel('Lookback Time (Myr)')
plt.ylabel('Mean Galactic Azimuth $\phi$ (rad)')
plt.legend()
plt.grid(True, linestyle=':', alpha=0.6)
plt.tight_layout()
plt.show()