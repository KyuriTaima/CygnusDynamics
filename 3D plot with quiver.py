import numpy as np
import pandas as pd
import plotly.graph_objects as go

# Data import
objects_data = pd.read_csv('Cygnus_Objects_Datas_Schönrich_239_8.3.csv')
objects_data = objects_data.dropna(subset=['Distance_pc', 'U_pec_kms', 'V_pec_kms', 'W_pec_kms'])
object_names = objects_data['Object_Name'].values
distances_pc = objects_data['Distance_pc'].values
l_deg = objects_data['l_deg'].values
b_deg = objects_data['b_deg'].values

u_pec = objects_data['U_pec_kms'].values
v_pec = objects_data['V_pec_kms'].values
w_pec = objects_data['W_pec_kms'].values

# Paramètres Galactiques
R0 = 8300 
l_rad, b_rad = np.radians(l_deg), np.radians(b_deg)

# CALCUL DES COORDONNÉES 3D (Centrées sur OB2)
x_gal = R0 - distances_pc * np.cos(l_rad) * np.cos(b_rad)
y_gal = distances_pc * np.sin(l_rad) * np.cos(b_rad)
z_gal = distances_pc * np.sin(b_rad)

# On trouve l'index de Group E
index_origin = list(object_names).index("Group E") 
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

fig.show()