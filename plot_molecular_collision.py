import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap

# --- Figure Initialization ---
fig, ax = plt.subplots(figsize=(10, 8), facecolor='#0B0C10')
ax.set_facecolor('#0B0C10')

# Define a professional scientific color palette
space_black = '#0B0C10'
halo_blue = '#1F2833'
filament_gray = '#C5C6C7'
accent_cyan = '#66FCF1'
text_white = '#FFFFFF'
flow_red = '#FF4F4F'

# --- 1. Drawing the Ambient Halo (Right Side) ---
# Create a smooth gradient interface for the gas halo boundary
gradient = np.linspace(0, 1, 100).reshape(1, -1)
ax.imshow(gradient, cmap=LinearSegmentedColormap.from_list('halo', [space_black, '#1A2536']), 
          extent=[4, 4.5, 0, 10], aspect='auto', zorder=1)

# Solid background for the rest of the dense halo
halo_rect = patches.Rectangle((4.5, 0), 5.5, 10, linewidth=0, facecolor='#1A2536', zorder=1)
ax.add_patch(halo_rect)

# --- 2. Drawing the Dense Filament (Left Side) ---
# Main body of the cylinder
cylinder_body = patches.Rectangle((1, 3.5), 3, 3, linewidth=1.5, edgecolor=filament_gray, facecolor='#2B2D31', zorder=3)
ax.add_patch(cylinder_body)

# Back ellipse of the cylinder
ellipse_back = patches.Ellipse((1, 5), 0.6, 3, linewidth=1.5, edgecolor=filament_gray, facecolor='#202225', zorder=2)
ax.add_patch(ellipse_back)

# Front ellipse of the cylinder (the nose entering the interface)
ellipse_front = patches.Ellipse((4, 5), 0.6, 3, linewidth=1.5, edgecolor=filament_gray, facecolor='#35383E', zorder=4)
ax.add_patch(ellipse_front)

# --- 3. Interaction Boundary & Vectors ---
# CORRECTION: Thickened Magnetic Field line (the main vertical boundary)
ax.axvline(x=4.0, color=accent_cyan, linestyle='-', linewidth=3.0, alpha=0.8, zorder=2)

# Velocity / Motion Vector of the filament
ax.annotate('', xy=(4.8, 5.0), xytext=(3.4, 5.0),
            arrowprops=dict(arrowstyle="->", color=accent_cyan, lw=2.5, mutation_scale=20), zorder=5)

# CORRECTION: Both arrows now represent the Particle Flow pointing towards the cylinder
# Top compressing particle flow
ax.annotate('', xy=(3.9, 6.7), xytext=(3.9, 8.5), arrowprops=dict(arrowstyle="->", color=flow_red, lw=1.8, mutation_scale=15))
# Bottom compressing particle flow
ax.annotate('', xy=(3.9, 3.3), xytext=(3.9, 1.5), arrowprops=dict(arrowstyle="->", color=flow_red, lw=1.8, mutation_scale=15))

# --- 4. Annotations & Scientific Labels ---
# Densities
ax.text(2.2, 5.0, r'$n = 10^5\ \mathrm{cm}^{-3}$', color=text_white, fontsize=12, fontweight='bold', ha='center', va='center', zorder=5)
ax.text(6.5, 8.5, r'$n = 50\ \mathrm{cm}^{-3}$', color=text_white, fontsize=14, fontweight='bold', ha='center')

# Scales and Dimensions
# Filament diameter (0.3 pc)
ax.annotate('', xy=(0.5, 3.5), xytext=(0.5, 6.5), arrowprops=dict(arrowstyle="<->", color=filament_gray, lw=1.2))
ax.text(0.3, 5.0, r'$0.3\ \mathrm{pc}$', color=filament_gray, fontsize=11, va='center', ha='right')

# Halo width (l ?)
ax.annotate('', xy=(4.0, 0.8), xytext=(9.5, 0.8), arrowprops=dict(arrowstyle="<->", color=text_white, lw=1.2))
ax.text(6.75, 1.2, r'$l = ?$', color=text_white, fontsize=13, fontweight='bold', ha='center')

# CORRECTION: Labels alignment and assignment
# Magnetic Field label at the top of the main line
ax.text(4.0, 9.3, 'Magnetic Field', color=accent_cyan, fontsize=11, ha='center', va='bottom', fontweight='bold')

# Symmetrical Particle Flow labels for both compressive arrows
ax.text(3.7, 7.6, 'Particle\nflow', color=flow_red, fontsize=10, ha='right', va='center', fontweight='semibold')
ax.text(3.7, 2.4, 'Particle\nflow', color=flow_red, fontsize=10, ha='right', va='center', fontweight='semibold')

# Relative velocity
ax.text(4.4, 5.3, r'$v_{\mathrm{rel}} \approx 15\ \mathrm{km/s}$', color=accent_cyan, fontsize=11, fontweight='semibold')

# --- 5. Axes Adjustments ---
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)
ax.axis('off')  # Hide standard plot axes for a clean look

# Save the professional figure
plt.tight_layout()
output_filename = "filament_interaction_schema.png"
plt.savefig(output_filename, dpi=300, facecolor=fig.get_facecolor(), edgecolor='none')
plt.show()