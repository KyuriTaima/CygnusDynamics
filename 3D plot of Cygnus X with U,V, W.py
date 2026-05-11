import numpy as np
import pandas as pd
import plotly.graph_objects as go

# 1. IMPORTATION ET NETTOYAGE DES DONNÉES
objects_data = pd.read_csv('Cygnus_Objects_Datas_Schönrich_239_8.3.csv')

# SÉCURITÉ : On supprime les lignes où il manque des valeurs de distance ou de vitesse
# (Sinon Plotly 3D plante en essayant de dessiner le vide)
objects_data = objects_data.dropna(subset=['Distance_pc', 'U_pec_kms', 'V_pec_kms', 'W_pec_kms'])

# Extraction sécurisée en tableaux Numpy (avec .values)
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

# 2. CALCUL DES COORDONNÉES 3D (Centrées sur OB2)
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

# 3. PRÉPARATION DE LA FIGURE
marker_colors = ['blue' if 'Group' in str(name) else 'black' for name in object_names]
fig = go.Figure()

# Ajout des points spatiaux (Nuages et Amas)
fig.add_trace(go.Scatter3d(
    x=x_local, y=y_local, z=z_local,
    mode='markers+text', text=object_names, textfont=dict(size=8, color=marker_colors),
    marker=dict(size=4, color=marker_colors, opacity=0.9),
    name='Cygnus X Objects'
))

# AJOUT DES VECTEURS VITESSE (U, V, W)
velocity_scale = 10.0

fig.add_trace(go.Cone(
    x=x_local, y=y_local, z=z_local,
    u=v_pec * velocity_scale,  # Axe X local (Rotation)
    v=u_pec * velocity_scale,  # Axe Y local (Centre Galactique)
    w=w_pec * velocity_scale,  # Axe Z local (Hauteur)
    colorscale='Viridis',
    showscale=True,
    colorbar=dict(title="Velocity<br>Magnitude", x=0.85),
    name='Velocity Vectors',
    hoverinfo='text',
    hovertext=[f"{name}<br>U: {u:.1f} km/s<br>V: {v:.1f} km/s<br>W: {w:.1f} km/s" 
               for name, u, v, w in zip(object_names, u_pec, v_pec, w_pec)]
))

# 5. CONFIGURATION VISUELLE
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
fig.write_html("Carte_Cygnus_3D_cones.html")