import numpy as np
import pandas as pd
import plotly.graph_objects as go
import matplotlib.pyplot as plt

# Data import
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

# fig.show()

# Create an html file to save for interactive plot
fig.write_html("Carte_Cygnus_3D_quiver.html")

# Now we go back in time using the velocity vectors to see the past positions of the objects, we can use the same time_array_myr for this purpose
past_positions = []
target_times = [4.5, 8.5, 10, 15]  # Time steps in Myr for which we want to plot the past positions
for t in target_times:
    x_past = x_local - v_pec * t * velocity_to_pc_myr  # V affects X
    y_past = y_local - u_pec * t * velocity_to_pc_myr  # U affects Y
    z_past = z_local - w_pec * t * velocity_to_pc_myr  # W affects Z
    past_positions.append((x_past, y_past, z_past))

    # plot the past positions on a new figure each time for specific time steps: 4.5 Myr, 8.5 Myr, 10 Myr, 15 Myr
    plt.figure(figsize=(10, 8))
    # Display names of the objects at their past positions
    plt.scatter(x_past, y_past, color='blue', label=f'-{t:.1f} Myr', alpha=0.7)
    for j in range(len(object_names)):
        plt.text(x_past[j], y_past[j], object_names[j], fontsize=10, color='black', weight='bold')
    plt.xlabel('X (pc)')
    plt.ylabel('Y (pc)')
    plt.title(f'Positions of Cygnus X Objects in the Galactic plane at -{t:.1f} Myr')
    plt.legend()
    plt.grid(True, linestyle=':', alpha=0.6)
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
        hoverinfo='skip' # On évite que la souris s'accroche au plan
    ))
    fig.update_layout(
        title=f"3D Kinematics of Cygnus X for positions at -{t:.1f} Myr",
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
    fig.write_html(f"Carte_Cygnus_3D_quiver_{t:.1f}Myr.html")

