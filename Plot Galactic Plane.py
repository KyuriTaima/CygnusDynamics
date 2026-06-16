import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Import csv data
objects_data = pd.read_csv('Cygnus_Objects_Datas_Uncertainties.csv')

# Extract relevant columns
names = objects_data['Object_Name'].values
U_pec_kms = objects_data['U_pec_kms'].values
V_pec_kms = objects_data['V_pec_kms'].values
distances = objects_data['Distance_pc'].values
l_deg = objects_data['l_deg'].values
b_deg = objects_data['b_deg'].values

R0 = 8277 # pc

# Galactic cartesian coordinates
l_rad = np.radians(l_deg)
b_rad = np.radians(b_deg)
X_gal = R0 - distances * np.cos(l_rad) * np.cos(b_rad)
Y_gal = distances * np.sin(l_rad) * np.cos(b_rad)

# We set the origin of the local frame on W75N
index_W75N = list(names).index("W75N")  # Index of W75N in the names list
X_W75N = X_gal[index_W75N]
Y_W75N = Y_gal[index_W75N]

# Global angle of Group E in the galactic plane relative to the galactic center
beta_W75N = np.arctan2(Y_W75N, X_W75N)

# Frame transformation from Galactic to Local (centered on W75N) for the kinematic vectors
# Translation (Origin = W75N)
delta_X = X_gal - X_W75N
delta_Y = Y_gal - Y_W75N

# Rotation of the axis
# Y_local towards the galactic center (aligned with U_pec)
Y_local = -(delta_X * np.cos(beta_W75N) + delta_Y * np.sin(beta_W75N))

# X_local axis : Points in the direction of the galactic rotation (aligned with V_pec)
X_local = -delta_X * np.sin(beta_W75N) + delta_Y * np.cos(beta_W75N)


# Graph creation
# set arrow colors: Red for groups A, D, F; Green for groups B, C, E; blue for the rest

arrow_colors = ['green' if 'Group B' in str(nom) or 'Group E' in str(nom) or 'Group C' in str(nom) else 'red' if 'Group A' in str(nom) or 'Group D' in str(nom) or 'Group F' in str(nom) else 'blue' for nom in names] 

plt.figure(figsize=(10, 10))

# Velocity vectors
Q = plt.quiver(X_local, Y_local, V_pec_kms, U_pec_kms, 
           angles='xy', scale_units='xy', scale=0.2, 
           color=arrow_colors, width=0.003, headwidth=3.5, headlength=4, zorder=3)

plt.quiverkey(Q, X=0.70, Y=0.97, U=10, 
              label='10 km/s', labelpos='E', 
              coordinates='axes', fontproperties={'weight': 'bold', 'size': 10}, 
              color='black')

# Set the origin (W75N) with a distinct marker
plt.scatter(0, 0, marker='o', s=20, color='black', edgecolor='black', zorder=2, label="W75N (Origin)")

# Add Sun direction
# The sun is located at (R0,0) in the galactic frame
X_gal_sun = R0
Y_gal_sun = 0

# Translation to the local frame centered on W75N
delta_X_sun = X_gal_sun - X_W75N
delta_Y_sun = Y_gal_sun - Y_W75N

# Frame rotation to get the sun's position in the local frame
Y_local_sun = -(delta_X_sun * np.cos(beta_W75N) + delta_Y_sun * np.sin(beta_W75N))
X_local_sun = -delta_X_sun * np.sin(beta_W75N) + delta_Y_sun * np.cos(beta_W75N)

# Trace the dashed line from the origin to the sun
# The sun is at ~1600 pc, so the line will extend well beyond the frame
plt.plot([0, X_local_sun], [0, Y_local_sun], color='orange', linestyle='--', linewidth=1, zorder=1)

# Add object names
for i in range(len(names)):
    couleur_texte = 'darkred' if 'Group' in str(names[i]) else 'teal'
    # Adjust the text position for IRAS to avoid overlap
    if names[i] == "IRAS20290+4052":
        xytext = (-50, 6)
    elif names[i] == "Group A":
        xytext = (-10, 6)
    else:
        xytext = (5, 5)
# annotate if the object is not W75N or Group E or Group D
    if names[i] not in ["W75N", "Group E", "Group D"]:
        plt.annotate(
            names[i], 
            (X_local[i], Y_local[i]),      
            textcoords="offset points",    
            xytext=xytext,                 
            ha='left',                     
            fontsize=10,
            color=couleur_texte
        )
    else:
        plt.annotate(
            names[i], 
            (X_local[i], Y_local[i]),      
            textcoords="offset points",    
            xytext=xytext,                 
            ha='left',                     
            fontsize=10,
            fontweight='bold',
            color=couleur_texte
        )
    if names[i] != "W75N":
        # Add a small marker for object, except for W75N which is already marked
        plt.scatter(X_local[i], Y_local[i], marker='o', s=20, color='white',edgecolor='black', zorder=2)



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
    xy=(600, 0),           
    xytext=(500, 0),       
    arrowprops=dict(facecolor='black', width=2, headwidth=5, alpha=0.6, shrink=0.05),
    fontsize=10, fontweight='bold', color='black',
    ha='right', va='center', zorder=3,
    bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1)
)

# Set markers for the axes
plt.axhline(0, color='black', linewidth=0.8, linestyle='--')
plt.axvline(0, color='black', linewidth=0.8, linestyle='--')

plt.ylabel('X-X(W75N): Distance along Galactic Rotation (pc)')
plt.xlabel('Y-Y(W75N): Distance towards Galactic Center (pc)')
plt.title('Local Kinematics of Cygnus X (Centered on W75N)')

plt.grid(True, linestyle=':', alpha=0.7)
plt.legend(loc='upper right')

# Set equal ratio for x and y axes
plt.gca().set_aspect('equal', adjustable='box')
plt.xlim(-100, 700)
plt.ylim(400, -100)

plt.show()

# Add a figure showing the positions in the X,Z plane (height above/below the galactic plane)
Z_gal = objects_data['Z_gc_kpc'].values
W_pec_kms = objects_data['W_pec_kms'].values

plt.figure(figsize=(10, 10))

# Velocity vectors
Q = plt.quiver(X_local, Z_gal, V_pec_kms, W_pec_kms, 
           angles='xz', scale_units='xz', scale=0.2, 
           color=arrow_colors, width=0.003, headwidth=3.5, headlength=4, zorder=3)

plt.quiverkey(Q, X=0.70, Y=0.97, U=10, 
              label='10 km/s', labelpos='E', 
              coordinates='axes', fontproperties={'weight': 'bold', 'size': 10}, 
              color='black')

plt.scatter(0, Z_gal[7], marker='o', s=20, color='black', edgecolor='black', zorder=2, label="W75N")

for i in range(len(names)):
    couleur_texte = 'darkred' if 'Group' in str(names[i]) else 'teal'
    # Adjust the text position for IRAS to avoid overlap
    if names[i] == "IRAS20290+4052":
        xytext = (-50, 6)
    elif names[i] == "Group A":
        xytext = (-10, 6)
    else:
        xytext = (5, 5)
# annotate if the object is not W75N or Group E or Group D
    if names[i] not in ["W75N", "Group E", "Group D"]:
        plt.annotate(
            names[i], 
            (X_local[i], Z_gal[i]),      
            textcoords="offset points",    
            xytext=xytext,                 
            ha='left',                     
            fontsize=10,
            color=couleur_texte
        )
    else:
        plt.annotate(
            names[i], 
            (X_local[i], Z_gal[i]),      
            textcoords="offset points",    
            xytext=xytext,                 
            ha='left',                     
            fontsize=10,
            fontweight='bold',
            color=couleur_texte
        )
    if names[i] != "W75N":
        # Add a small marker for object, except for W75N which is already marked
        plt.scatter(X_local[i], Z_gal[i], marker='o', s=20, color='white',edgecolor='black', zorder=2)

plt.annotate(
    'To Galactic Center', 
    xy=(0, 20),             
    xytext=(0, 300),         
    arrowprops=dict(facecolor='black', width=2, headwidth=8, alpha=0.6, shrink=0.05),
    fontsize=10, fontweight='bold', color='black',
    ha='center', va='top', zorder=4,
    bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1)
)

# Set markers for the axes
plt.axhline(0, color='black', linewidth=0.8, linestyle='--')
plt.axvline(0, color='black', linewidth=0.8, linestyle='--')

plt.xlabel('X-X(W75N): Distance along Galactic Rotation (pc)')
plt.ylabel('Z: Distance towards Galactic North Pole (pc)')
plt.title('Local Kinematics of Cygnus X in the X-Z Plane (Centered on W75N)')

plt.grid(True, linestyle=':', alpha=0.7)
plt.legend(loc='upper right')

# Set equal ratio for x and y axes
plt.gca().set_aspect('equal', adjustable='box')
plt.xlim(-100, 700)
plt.ylim(400, -100)

plt.show()

