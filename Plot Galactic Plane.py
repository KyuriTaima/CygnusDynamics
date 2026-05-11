import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Import csv data
objects_data = pd.read_csv('Cygnus_Objects_Datas_Reid2019_236_8.15.csv')

# Extract relevant columns
names = objects_data['Object_Name'].values
U_pec_kms = objects_data['U_pec_kms'].values
V_pec_kms = objects_data['V_pec_kms'].values
distances = objects_data['Distance_pc'].values
l_deg = objects_data['l_deg'].values
b_deg = objects_data['b_deg'].values

R0 = 8150 # pc

# Galactic cartesian coordinates
l_rad = np.radians(l_deg)
b_rad = np.radians(b_deg)
X_gal = R0 - distances * np.cos(l_rad) * np.cos(b_rad)
Y_gal = distances * np.sin(l_rad) * np.cos(b_rad)

# We set the origin of the local frame on Cygnus OB 2 (Group E)
index_E = list(names).index("Group E") 
X_E = X_gal[index_E]
Y_E = Y_gal[index_E]

# Global angle of Group E in the galactic plane relative to the galactic center
beta_E = np.arctan2(Y_E, X_E)

# Frame transformation from Galactic to Local (centered on Group E) for the kinematic vectors
# Translation (Origin = Group E)
delta_X = X_gal - X_E
delta_Y = Y_gal - Y_E

# Rotation of the axis
# Y_local towards the galactic center (aligned with U_pec)
Y_local = -(delta_X * np.cos(beta_E) + delta_Y * np.sin(beta_E))

# X_local axis : Points in the direction of the galactic rotation (aligned with V_pec)
X_local = -delta_X * np.sin(beta_E) + delta_Y * np.cos(beta_E)


# Graph creation
arrow_colors = ['red' if 'Group' in str(nom) else 'cyan' for nom in names]

plt.figure(figsize=(10, 10))

# Velocity vectors
plt.quiver(X_local, Y_local, V_pec_kms, U_pec_kms, 
           angles='xy', scale_units='xy', scale=0.2, 
           color=arrow_colors, width=0.003, headwidth=3.5, headlength=4, zorder=3)

# Set the origin (Cygnus OB 2 / Group E) with a distinct marker
plt.scatter(0, 0, marker='o', s=50, color='black', edgecolor='black', zorder=2, label="Cygnus OB 2 (Origin)")

# Add Sun direction
# The sun is located at (R0,0) in the galactic frame
X_gal_sun = R0
Y_gal_sun = 0

# Translation to the local frame centered on Group E
delta_X_sun = X_gal_sun - X_E
delta_Y_sun = Y_gal_sun - Y_E

# Frame rotation to get the sun's position in the local frame
Y_local_sun = -(delta_X_sun * np.cos(beta_E) + delta_Y_sun * np.sin(beta_E))
X_local_sun = -delta_X_sun * np.sin(beta_E) + delta_Y_sun * np.cos(beta_E)

# Trace the dashed line from the origin to the sun
# The sun is at ~1600 pc, so the line will extend well beyond the frame
plt.plot([0, X_local_sun], [0, Y_local_sun], color='orange', linestyle='--', linewidth=1, zorder=1)

# Add object names
for i in range(len(names)):
    couleur_texte = 'darkred' if 'Group' in str(names[i]) else 'teal'
    plt.annotate(
        names[i], 
        (X_local[i], Y_local[i]),      
        textcoords="offset points",    
        xytext=(6, 6),                 
        ha='left',                     
        fontsize=10,
        fontweight='bold',
        color=couleur_texte
    )


# Add galactic indicators

# Arrow towards the galactic center
plt.annotate(
    'To Galactic Center', 
    xy=(0, 380),             
    xytext=(0, 300),         
    arrowprops=dict(facecolor='black', width=2, headwidth=8, alpha=0.6, shrink=0.05),
    fontsize=10, fontweight='bold', color='black',
    ha='center', va='top', zorder=4,
    bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1)
)

# Arrow in the direction of the galactic rotation
plt.annotate(
    'Galactic Rotation', 
    xy=(380, 0),           
    xytext=(300, 0),       
    arrowprops=dict(facecolor='black', width=2, headwidth=5, alpha=0.6, shrink=0.05),
    fontsize=10, fontweight='bold', color='black',
    ha='right', va='center', zorder=3,
    bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1)
)

# Set markers for the axes
plt.axhline(0, color='black', linewidth=0.8, linestyle='--')
plt.axvline(0, color='black', linewidth=0.8, linestyle='--')

plt.xlabel('Distance along Galactic Rotation (pc)')
plt.ylabel('Distance towards Galactic Center (pc)')
plt.title('Local Kinematics of Cygnus X (Centered on Cygnus OB2 / Group E)')

plt.grid(True, linestyle=':', alpha=0.7)
plt.legend(loc='upper right')

# Set equal ratio for x and y axes
plt.gca().set_aspect('equal', adjustable='box')
plt.xlim(-400, 400)
plt.ylim(-400, 400)

plt.show()