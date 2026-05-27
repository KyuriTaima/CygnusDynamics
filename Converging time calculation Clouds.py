# Calculation of the converging time of the molecular clouds
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import combinations

# --- Data Loading and Filtering ---
data_file = 'Cygnus_Objects_Datas_Uncertainties.csv'
objects_data = pd.read_csv(data_file)

# Filter data to isolate the BEC supergroup
target_clusters = [ "W75N", "DR21", "DR20", "IRAS20290+4052"]
bec_data = objects_data[objects_data['Object_Name'].isin(target_clusters)].reset_index(drop=True)
cluster_names = bec_data['Object_Name'].values

# --- Parameter Extraction and Unit Conversion ---
# Convert Galactocentric spatial coordinates from kpc to pc for precision
X_0 = bec_data['X_gc_pc'].values * 1000  
Y_0 = bec_data['Y_gc_pc'].values * 1000
Z_0 = bec_data['Z_gc_pc'].values * 1000

# Galactocentric velocity components in km/s
Vx = bec_data['vx_gc_kms'].values
Vy = bec_data['vy_gc_kms'].values
Vz = bec_data['vz_gc_kms'].values

# --- Kinematic Traceback Setup ---
# Time array from 0 to 20 Myr with a 0.1 Myr resolution
time_array_myr = np.arange(0.0, 20.1, 0.1)
mean_distances_array = []

# Generate all unique pairs of clusters to compute relative distances
cluster_pairs = list(combinations(range(len(cluster_names)), 2))

# Conversion factor: 1 km/s equals approximately 1.022 pc/Myr
velocity_to_pc_myr = 1.022

# --- Traceback Execution ---
for t in time_array_myr:
    # Compute past positions assuming linear ballistic motion
    X_t = X_0 + (Vx * t * velocity_to_pc_myr)
    Y_t = Y_0 + (Vy * t * velocity_to_pc_myr)
    Z_t = Z_0 + (Vz * t * velocity_to_pc_myr)
    
    # Compute 3D Euclidean distance for each pair at epoch 't'
    distances_at_epoch = []
    for i, j in cluster_pairs:
        dist = np.sqrt((X_t[i] - X_t[j])**2 + (Y_t[i] - Y_t[j])**2 + (Z_t[i] - Z_t[j])**2)
        distances_at_epoch.append(dist)
        
    # Store the mean relative distance of the supergroup at this epoch
    mean_distances_array.append(np.mean(distances_at_epoch))

mean_distances_array = np.array(mean_distances_array)

# --- Kinematic Age Determination (Time of Closest Approach) ---
min_distance_index = np.argmin(mean_distances_array)
kinematic_age = time_array_myr[min_distance_index]
minimum_mean_distance = mean_distances_array[min_distance_index]

# Output results to console
print(f"--- Kinematic Prevision Results for {', '.join(cluster_names)} ---")
print(f"Time of Closest Approach (TCA) : {kinematic_age:.1f} Myr")
print(f"Minimum Mean Relative Distance : {minimum_mean_distance:.1f} pc")
print("-" * 50)

# --- Visualization ---
plt.figure(figsize=(9, 6))
plt.plot(time_array_myr, mean_distances_array, lw=2.5, color='royalblue', label='Mean Relative Distance')

