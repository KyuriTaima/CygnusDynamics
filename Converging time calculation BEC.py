import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import combinations

# Data loading
data_file = 'Cygnus_Objects_Datas_Uncertainties.csv'
objects_data = pd.read_csv(data_file)

# Filter data to isolate the BEC supergroup
target_clusters = ["Group B", "Group E", "Group C"]
bec_data = objects_data[objects_data['Object_Name'].isin(target_clusters)].reset_index(drop=True)
cluster_names = bec_data['Object_Name'].values

# Convert Galactocentric spatial coordinates from kpc to pc for precision
X_0_nom = bec_data['X_gc_kpc'].values * 1000  
Y_0_nom = bec_data['Y_gc_kpc'].values * 1000
Z_0_nom = bec_data['Z_gc_kpc'].values * 1000

# Galactocentric velocity components in km/s
Vx_nom = bec_data['vx_gc_kms'].values
Vy_nom = bec_data['vy_gc_kms'].values
Vz_nom = bec_data['vz_gc_kms'].values

# Uncertainties
e_X_0 = bec_data['X_gc_err_kpc'].values * 1000
e_Y_0 = bec_data['Y_gc_err_kpc'].values * 1000
e_Z_0 = bec_data['Z_gc_err_kpc'].values * 1000
e_Vx = bec_data['vx_gc_err_kms'].values
e_Vy = bec_data['vy_gc_err_kms'].values
e_Vz = bec_data['vz_gc_err_kms'].values

# Kinematic traceback parameters
time_array_myr = np.arange(0.0, 20.1, 0.1)
cluster_pairs = list(combinations(range(len(cluster_names)), 2))
velocity_to_pc_myr = 1.022

# Monte Carlo Simulation for Uncertainty Propagation
N_simulations = 10000  # Number of Monte Carlo simulations
all_mean_distances = np.zeros((N_simulations, len(time_array_myr)))
all_tcas = []
all_min_distances = []

print(f"Lancement de {N_simulations} simulations de Monte-Carlo...")

for sim in range(N_simulations):
    # Normal random sampling for initial positions and velocities based on uncertainties
    X_0 = np.random.normal(X_0_nom, e_X_0)
    Y_0 = np.random.normal(Y_0_nom, e_Y_0)
    Z_0 = np.random.normal(Z_0_nom, e_Z_0)
    Vx = np.random.normal(Vx_nom, e_Vx)
    Vy = np.random.normal(Vy_nom, e_Vy)
    Vz = np.random.normal(Vz_nom, e_Vz)
    
    # Execute the kinematic traceback for each time step
    sim_distances = []
    for t in time_array_myr:
        X_t = X_0 - (Vx * t * velocity_to_pc_myr)
        Y_t = Y_0 - (Vy * t * velocity_to_pc_myr)
        Z_t = Z_0 - (Vz * t * velocity_to_pc_myr)
        
        # Calculate mean relative distances between all cluster pairs at this time step
        dists = []
        for i, j in cluster_pairs:
            dist = np.sqrt((X_t[i] - X_t[j])**2 + (Y_t[i] - Y_t[j])**2 + (Z_t[i] - Z_t[j])**2)
            dists.append(dist)
        sim_distances.append(np.mean(dists))
        
    all_mean_distances[sim, :] = sim_distances
    
    # Get the time of closest approach (TCA) and the corresponding minimum mean distance
    min_idx = np.argmin(sim_distances)
    all_tcas.append(time_array_myr[min_idx])
    all_min_distances.append(sim_distances[min_idx])

# Statistical analysis of the Monte Carlo results
# Distances
median_distances = np.median(all_mean_distances, axis=0)
lower_bound_dist = np.percentile(all_mean_distances, 16, axis=0) # -1 sigma
upper_bound_dist = np.percentile(all_mean_distances, 84, axis=0) # +1 sigma

# TCA and minimum mean distance statistics
tca_median = np.median(all_tcas)
tca_std = np.std(all_tcas)
min_dist_median = np.median(all_min_distances)
min_dist_std = np.std(all_min_distances)

# Output results to console
print(f"\n--- Kinematic Traceback Results for {', '.join(cluster_names)} ---")
print(f"Time of Closest Approach (TCA) : -{tca_median:.1f} ± {tca_std:.1f} Myr")
print(f"Minimum Mean Relative Distance : {min_dist_median:.1f} ± {min_dist_std:.1f} pc")
print("-" * 50)

# Vizualisation of the results
plt.figure(figsize=(10, 6))

# Add y axis label
plt.ylabel('Mean Relative Distance (pc)', fontsize=12, fontweight='bold')
# Add x axis label
plt.xlabel('Time in the past (Myr)', fontsize=12, fontweight='bold')
# Add title
plt.title('Kinematic Traceback of the E-BC Supergroup with Monte-Carlo Uncertainties', fontsize=14, fontweight='bold')

# Draw the uncertainty band (1 sigma) around the median trajectory
plt.fill_between(time_array_myr, lower_bound_dist, upper_bound_dist, color='royalblue', alpha=0.3, label=r'1$\sigma$ Uncertainty Band')

# Draw the median trajectory of the mean relative distance over time
plt.plot(time_array_myr, median_distances, lw=2.5, color='royalblue', label='Median Relative Distance')

# Add a vertical dashed line at the time of closest approach with uncertainty
plt.axvline(x=tca_median, color='red', linestyle='--', lw=2, label=f'TCA: -{tca_median:.1f} $\pm$ {tca_std:.1f} Myr')
plt.axvspan(tca_median - tca_std, tca_median + tca_std, color='red', alpha=0.15) # Zone d'incertitude du TCA

# Add minimal mean relative distance point
plt.scatter(tca_median, min_dist_median, color='red', s=100, zorder=5, label=f'Min Dist: {min_dist_median:.1f} $\pm$ {min_dist_std:.1f} pc')  

plt.grid(True, linestyle=':', alpha=0.6)
plt.legend()
plt.tight_layout()
plt.xlim(0.0, 12.5) # Limit the x-axis to 12.5 Myr for better visualization of the TCA region
plt.show()