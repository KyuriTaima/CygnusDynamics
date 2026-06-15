# Linear traceback for absolute values

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import matplotlib.pyplot as plt

# Data import
theta_0 = 234.8
objects_data = pd.read_csv('Cygnus_Objects_Datas_Uncertainties.csv')
objects_data = objects_data.dropna(subset=['Distance_pc', 'U_pec_kms', 'V_pec_kms', 'W_pec_kms'])
object_names = objects_data['Object_Name'].values
distances_pc = objects_data['Distance_pc'].values
l_deg = objects_data['l_deg'].values
b_deg = objects_data['b_deg'].values

u_pec = objects_data['U_pec_kms'].values
v_pec = objects_data['V_pec_kms'].values
w_pec = objects_data['W_pec_kms'].values

# Création d'un pas de temps de 20 Myr avec une résolution de 0.1 Myr
time_array_myr = np.arange(0.0, 20.1, 0.1)

# Conversion factor: 1 km/s equals approximately 1.022 pc/Myr
velocity_to_pc_myr = 1.022

# Paramètres Galactiques
R0 = 8277 
l_rad, b_rad = np.radians(l_deg), np.radians(b_deg)

# CALCUL DES COORDONNÉES 3D (Centrées sur W75N)
x_gal = R0 - distances_pc * np.cos(l_rad) * np.cos(b_rad)
y_gal = distances_pc * np.sin(l_rad) * np.cos(b_rad)
z_gal = distances_pc * np.sin(b_rad)

# On trouve l'index de Group E
index_origin = list(object_names).index("W75N") 
x_origin, y_origin, z_origin = x_gal[index_origin], y_gal[index_origin], z_gal[index_origin]
beta_angle = np.arctan2(y_origin, x_origin)

# Translation puis Rotation
dx, dy, dz = x_gal - x_origin, y_gal - y_origin, z_gal - z_origin
y_local = -(dx * np.cos(beta_angle) + dy * np.sin(beta_angle))
x_local = -dx * np.sin(beta_angle) + dy * np.cos(beta_angle)
z_local = dz

# PRÉPARATION DE LA FIGURE
marker_colors = ['blue' if 'Group' in str(name) else 'black' for name in object_names]
fig = go.Figure()

# Ajout des points spatiaux (Nuages et Amas)
fig.add_trace(go.Scatter3d(
    x=x_local, y=y_local, z=z_local,
    mode='markers+text', text=object_names, textfont=dict(size=8, color=marker_colors),
    marker=dict(size=4, color=marker_colors, opacity=0.9),
    name='Cygnus X Objects'
))

# AJOUT DES VECTEURS VITESSE
velocity_scale = 5.0

# On boucle sur chaque objet pour tracer sa ligne de vitesse
for i in range(len(object_names)):
    # Calcul des coordonnées de la pointe de la flèche
    end_x = x_local[i] + v_pec[i] * velocity_scale  # V va sur X
    end_y = y_local[i] + u_pec[i] * velocity_scale  # U va sur Y
    end_z = z_local[i] + w_pec[i] * velocity_scale  # W va sur Z
    
    # Calcul de la norme pour l'affichage au survol
    v_norm = np.sqrt(u_pec[i]**2 + v_pec[i]**2 + w_pec[i]**2)
    
    fig.add_trace(go.Scatter3d(
        x=[x_local[i], end_x],
        y=[y_local[i], end_y],
        z=[z_local[i], end_z],
        mode='lines',
        line=dict(
            color='gold', 
            width=5       
        ),
        showlegend=False,
        hoverinfo='text',
        hovertext=f"{object_names[i]} Velocity<br>|V| = {v_norm:.1f} km/s<br>U:{u_pec[i]:.1f} V:{v_pec[i]:.1f} W:{w_pec[i]:.1f}"
    ))

# CONFIGURATION VISUELLE
fig.update_layout(
    title="3D Kinematics of Cygnus X (Positions + Velocity Vectors)",
    scene=dict(
        xaxis_title='Rotation (pc)', 
        yaxis_title='Center (pc)', 
        zaxis_title='Height Z (pc)',
        aspectmode='data', 
        bgcolor='rgb(10, 10, 20)'
    ),
    paper_bgcolor='rgb(10, 10, 20)', font=dict(color='white'),
    margin=dict(l=0, r=0, b=0, t=40)
)

# fig.show()

# Create an html file to save for interactive plot
fig.write_html("Carte_Cygnus_3D_quiver.html")

# ==============================================================================
# RÉGRESSION LINÉAIRE DU BRAS SPIRAL (Inclinaison 11.4°)
# ==============================================================================
pitch_angle_deg = 11.4
cot_i = 1 / np.tan(np.radians(pitch_angle_deg))
arm_width_pc = 300 # Largeur du bras

# Données d'impact des amas : [W75N, Groupe E, Groupe D]
t_impacts = np.array([0.0, -4.5, -9.8])
x_impacts = np.array([0.0, 399.0, 730.0])
y_impacts = np.array([0.0, 54.0, 187.0])

# On calcule la projection spatiale D
D_proj = x_impacts - (y_impacts * cot_i)

# Régression linéaire : D_proj = W * t + X0
# np.polyfit(x, y, 1) retourne [pente, ordonnée_origine] pour une fonction de degré 1
W_rel_pc_myr, X0_pc = np.polyfit(t_impacts, D_proj, 1)

print(f"--- RÉSULTATS DE LA RÉGRESSION ---")
print(f"Vitesse relative de l'onde (W) : {W_rel_pc_myr:.1f} pc/Myr")
print(f"Décalage spatial à t=0 (X0) : {X0_pc:.1f} pc")
print(f"----------------------------------")

past_positions = []
target_times = [0.0, 4.5, 9.8, 10, 15]  # Time steps in Myr

# ==============================================================================
# PARAMÈTRES DU BRAS SPIRAL (Modèle de Reid 2019 : i = 11.4°)
# ==============================================================================
pitch_angle_deg = -11.4
cot_i = 1 / np.tan(np.radians(pitch_angle_deg))
W_rel_pc_myr = 21.4  # Vitesse relative de l'onde calculée pour l'impact

past_positions = []
target_times = [0.0, 4.5, 9.8, 10, 15]  # Time steps in Myr


for t_past in target_times:
    t_phys = -t_past # Le temps recule (valeurs négatives)
    
    x_past = x_local + v_pec * t_phys * velocity_to_pc_myr
    y_past = y_local + u_pec * t_phys * velocity_to_pc_myr
    z_past = z_local + w_pec * t_phys * velocity_to_pc_myr
    past_positions.append((x_past, y_past, z_past))

    plt.figure(figsize=(10, 8))
    
    # Trace des amas
    marker_colors = ['red' if 'Group F' in str(name) or 'Group A' in str(name) or 'Group D' in str(name) else 'green' if 'Group B' in str(name) or 'Group E' in str(name) or 'Group C' in str(name) else 'blue' for name in object_names]
    plt.scatter(x_past, y_past, color=marker_colors, label=f'Amas à t={t_phys:.1f} Myr', alpha=0.7, zorder=5)
    for j in range(len(object_names)):
        plt.text(x_past[j], y_past[j] + 15, object_names[j], fontsize=10, color='black', weight='bold', zorder=5)
    
# -------------------------------------------------------------
    # TRACÉ DU BRAS SPIRAL OPTIMISÉ PAR RÉGRESSION
    # -------------------------------------------------------------
    y_arm_line = np.linspace(-300, 800, 100) 
    
    # Équation ajustée avec X0
    x_arm_center = y_arm_line * cot_i + (W_rel_pc_myr * t_phys) + X0_pc
    
    width_x_offset = (arm_width_pc / 2) / np.cos(np.radians(pitch_angle_deg))
    x_arm_left = x_arm_center - width_x_offset
    x_arm_right = x_arm_center + width_x_offset
    
    # Tracé de l'épaisseur et du centre
    plt.fill_betweenx(y_arm_line, x_arm_left, x_arm_right, color='magenta', alpha=0.15, label=f'Épaisseur Bras ({arm_width_pc} pc)')
    plt.plot(x_arm_center, y_arm_line, color='magenta', linestyle='--', linewidth=2.5, label=f'Centre Ajusté (W={W_rel_pc_myr:.1f} pc/Myr)')
    # -------------------------------------------------------------

    plt.xlabel('Distance along Galactic Rotation (pc)')
    plt.ylabel('Distance towards Galactic Center (pc)')
    plt.title(f'Positions of Cygnus X Objects and Spiral Arm at t={t_phys:.1f} Myr')
    plt.gca().set_aspect('equal', adjustable='box')
    
    # Axes avec Y inversé (Centre galactique vers le bas)
    plt.xlim(-100, 1000)
    plt.ylim(500, -100) 

    plt.annotate(
        'To Galactic Center', 
        xy=(40, 480),            
        xytext=(40, 400),        
        arrowprops=dict(facecolor='black', width=2, headwidth=8, alpha=0.6, shrink=0.05),
        fontsize=10, fontweight='bold', color='black',
        ha='center', va='top', zorder=4,
        bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1)
    )

    plt.annotate(
        'Galactic Rotation', 
        xy=(950, 450),           
        xytext=(850, 450),       
        arrowprops=dict(facecolor='black', width=2, headwidth=5, alpha=0.6, shrink=0.05),
        fontsize=10, fontweight='bold', color='black',
        ha='right', va='center', zorder=3,
        bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1)
    )
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.legend(loc='upper right')
    plt.show()

    # Also export a 3D version of the past positions for this time step
    marker_colors = ['blue' if 'Group' in str(name) else 'black' for name in object_names]
    fig = go.Figure()
    fig.add_trace(go.Scatter3d(
        x=x_past, y=y_past, z=z_past,
        mode='markers+text', text=object_names, textfont=dict(size=8, color=marker_colors),
        marker=dict(size=4, color=marker_colors, opacity=0.9),
        name='Cygnus X Objects'))
    plane_size = max(np.max(np.abs(x_past)), np.max(np.abs(y_past))) * 1.2
    
    fig.add_trace(go.Mesh3d(
        x=[-plane_size, plane_size, plane_size, -plane_size],
        y=[-plane_size, -plane_size, plane_size, plane_size],
        z=[0, 0, 0, 0],
        color='cyan',
        opacity=0.15,
        name='Galactic Plane (Z=0)',
        hoverinfo='skip' 
    ))
    fig.update_layout(
        title=f"3D Kinematics of Cygnus X for positions at -{t_past:.1f} Myr",
        scene=dict(
            xaxis_title='Rotation (pc)', 
            yaxis_title='Center (pc)', 
            zaxis_title='Height Z (pc)',
            aspectmode='data', 
            bgcolor='rgb(10, 10, 20)'
        ),
        paper_bgcolor='rgb(10, 10, 20)', font=dict(color='white'),
        margin=dict(l=0, r=0, b=0, t=40)
    )
    fig.write_html(f"Carte_Cygnus_3D_quiver_{-t_phys:.1f}Myr.html")