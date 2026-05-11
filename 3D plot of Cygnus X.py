import numpy as np
import pandas as pd
import plotly.graph_objects as go

# Import dataset and extract relevant columns
objects_data = pd.read_csv('Cygnus_Objects_Datas_Schönrich_239_8.3.csv')
object_names = objects_data['Object_Name'].values
distances_pc = objects_data['Distance_pc'].values
l_deg = objects_data['l_deg'].values
b_deg = objects_data['b_deg'].values

# Distance from the Sun to the Galactic Center in parsecs
R0 = 8300 

# Convert galactic coordinates from degrees to radians
l_rad = np.radians(l_deg)
b_rad = np.radians(b_deg)

# Calculate global Galactocentric coordinates (X, Y, Z)
x_gal = R0 - distances_pc * np.cos(l_rad) * np.cos(b_rad)
y_gal = distances_pc * np.sin(l_rad) * np.cos(b_rad)
z_gal = distances_pc * np.sin(b_rad)

# Define Cygnus OB2 (Group E) as the local origin
index_origin = list(object_names).index("Group E") 
x_origin = x_gal[index_origin]
y_origin = y_gal[index_origin]
z_origin = z_gal[index_origin]

# Calculate the rotation angle of the origin relative to the Galactic Center
beta_angle = np.arctan2(y_origin, x_origin)

# Apply translation to center the frame on Cygnus OB2
dx = x_gal - x_origin
dy = y_gal - y_origin
dz = z_gal - z_origin

# Apply rotation to align local axes (X: Rotation, Y: Galactic Center)
y_local = -(dx * np.cos(beta_angle) + dy * np.sin(beta_angle))
x_local = -dx * np.sin(beta_angle) + dy * np.cos(beta_angle)
z_local = dz

# Assign colors based on object type (red for clusters, cyan for clouds)
marker_colors = ['blue' if 'Group' in str(name) else 'black' for name in object_names]

# Initialize the 3D figure
fig = go.Figure()

# Add scatter plot for all objects
fig.add_trace(go.Scatter3d(
    x=x_local, 
    y=y_local, 
    z=z_local,
    mode='markers+text',
    text=object_names,
    marker=dict(size=5, color=marker_colors, opacity=0.8),
    name='Cygnus X Objects'
))

# Configure scene layout and axis titles
fig.update_layout(
    title="Internal 3D Structure of Cygnus X (Centered on OB2)",
    scene=dict(
        xaxis_title='Towards Galactic Rotation (pc)',
        yaxis_title='Towards Galactic Center (pc)',
        zaxis_title='Galactic Height Z (pc)',
        aspectmode='data' # Ensures true physical proportions in 3D space
    ),
    margin=dict(l=0, r=0, b=0, t=40)
)

# Display the interactive plot
fig.show()