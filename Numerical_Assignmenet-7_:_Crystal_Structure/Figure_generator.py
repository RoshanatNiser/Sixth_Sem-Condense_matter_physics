"""
generate_fig1.py
----------------
Generates Figure 1 for the 2D Crystallographic Point Groups assignment.
Produces:  Fig-1.png  (300 dpi, for Overleaf upload)
           Fig-1.pdf  (vector, also usable in Overleaf)

Usage:
    python generate_fig1.py

Dependencies: matplotlib, numpy (both standard in any scientific Python install)
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# ------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------

def setup(ax, label):
    """Clean axes with a border box and bold title."""
    ax.set_xlim(-1.3, 1.3)
    ax.set_ylim(-1.3, 1.3)
    ax.set_aspect('equal')
    ax.axis('off')
    border = mpatches.FancyBboxPatch(
        (-1.25, -1.25), 2.5, 2.5,
        boxstyle="square,pad=0",
        edgecolor='black', facecolor='white', linewidth=1.2
    )
    ax.add_patch(border)
    ax.set_title(label, fontsize=10, fontweight='bold', pad=3)


def dot(ax, x, y, color='black', s=55):
    """Draw a filled dot at (x, y)."""
    ax.scatter([x], [y], color=color, s=s, zorder=6)


def rect_motif(ax, w=1.0, h=0.62):
    """Draw a rectangle motif."""
    r = mpatches.Rectangle(
        (-w/2, -h/2), w, h,
        edgecolor='black', facecolor='#e8f4f8', lw=1.4
    )
    ax.add_patch(r)


def sq_motif(ax, s=0.9):
    """Draw a square motif."""
    r = mpatches.Rectangle(
        (-s/2, -s/2), s, s,
        edgecolor='black', facecolor='#e8f4f8', lw=1.4
    )
    ax.add_patch(r)


def tri_motif(ax, r=0.85):
    """Draw an equilateral triangle motif (point up)."""
    angs = [np.pi/2, np.pi/2 + 2*np.pi/3, np.pi/2 + 4*np.pi/3]
    xs = [r * np.cos(a) for a in angs]
    ys = [r * np.sin(a) for a in angs]
    ax.fill(xs, ys, edgecolor='black', facecolor='#e8f8e8', lw=1.4)


def hex_motif(ax, r=0.88):
    """Draw a regular hexagon motif."""
    angs = np.linspace(0, 2*np.pi, 7)
    ax.fill(
        r * np.cos(angs), r * np.sin(angs),
        edgecolor='black', facecolor='#f8f0e8', lw=1.4
    )


# ------------------------------------------------------------------
# Build figure
# ------------------------------------------------------------------

fig, axes = plt.subplots(2, 5, figsize=(14, 6))
axes = axes.flatten()

# ---- C1: single asymmetric dot ----
ax = axes[0]
setup(ax, r'$C_1$')
rect_motif(ax)
dot(ax, 0.32, 0.22)
ax.text(0, -1.15, 'oblique, 1 dot\n(no symm.)',
        ha='center', fontsize=7.5, style='italic')

# ---- C2: two dots related by 180-deg rotation ----
ax = axes[1]
setup(ax, r'$C_2$')
rect_motif(ax)
dot(ax,  0.35,  0.22)
dot(ax, -0.35, -0.22)
ax.annotate('', xy=(-0.30, -0.18), xytext=(0.30, 0.18),
            arrowprops=dict(
                arrowstyle='->,head_width=0.15',
                color='crimson', lw=1.0,
                connectionstyle='arc3,rad=0.5'
            ))
ax.text(0, -1.15, 'rect., 2 dots\n$C_2$ related',
        ha='center', fontsize=7.5, style='italic')

# ---- C3: chiral triangle (dots rotated off mirror lines) ----
ax = axes[2]
setup(ax, r'$C_3$')
tri_motif(ax)
r = 0.48
offset = 0.35   # chiral offset so dots are NOT on the altitudes
for i in range(3):
    ang = np.pi/2 + offset + i * 2*np.pi/3
    dot(ax, r * np.cos(ang), r * np.sin(ang))
ax.text(0, -1.15, 'hex., 3 dots\nchiral (no mirror)',
        ha='center', fontsize=7.5, style='italic')

# ---- C4: pinwheel (dots rotated off mirror lines) ----
ax = axes[3]
setup(ax, r'$C_4$')
sq_motif(ax)
r = 0.42
offset = 0.22   # rotate so dots land between the 4 mirror lines
for i in range(4):
    ang = np.pi/4 + offset + i * np.pi/2
    dot(ax, r * np.cos(ang), r * np.sin(ang))
ax.text(0, -1.15, 'sq., pinwheel\nno mirror',
        ha='center', fontsize=7.5, style='italic')

# ---- C6: pinwheel hexagon ----
ax = axes[4]
setup(ax, r'$C_6$')
hex_motif(ax)
r = 0.52
offset = 0.22
for i in range(6):
    ang = i * np.pi/3 + offset
    dot(ax, r * np.cos(ang), r * np.sin(ang))
ax.text(0, -1.15, 'hex., 6 dots\nchiral pinwheel',
        ha='center', fontsize=7.5, style='italic')

# ---- D1: one mirror, dot ON vertical mirror only ----
ax = axes[5]
setup(ax, r'$D_1$')
rect_motif(ax)
ax.axvline(0, color='red', lw=1.0, ls='--', alpha=0.7,
           ymin=0.08, ymax=0.92)
dot(ax, 0.0, 0.25)
ax.text(0.08, 0.55, r'$\sigma$', color='red', fontsize=9)
ax.text(0, -1.15, 'rect., dot on\none mirror only',
        ha='center', fontsize=7.5, style='italic')

# ---- D2: two mirrors, four symmetric dots ----
ax = axes[6]
setup(ax, r'$D_2$')
rect_motif(ax)
ax.axvline(0, color='red',  lw=1.0, ls='--', alpha=0.6,
           ymin=0.08, ymax=0.92)
ax.axhline(0, color='blue', lw=1.0, ls='--', alpha=0.6,
           xmin=0.08, xmax=0.92)
for x, y in [(0.33, 0.22), (-0.33, 0.22),
             (0.33, -0.22), (-0.33, -0.22)]:
    dot(ax, x, y)
ax.text(0, -1.15, 'rect., 4 dots\n2 mirrors + $C_2$',
        ha='center', fontsize=7.5, style='italic')

# ---- D3: dots ON the three altitudes ----
ax = axes[7]
setup(ax, r'$D_3$')
tri_motif(ax)
r_mirror = 0.82
r_dot    = 0.46
for i in range(3):
    ang = np.pi/2 + i * 2*np.pi/3
    # mirror line
    ax.plot([0, r_mirror * np.cos(ang)],
            [0, r_mirror * np.sin(ang)],
            'r--', lw=0.9, alpha=0.65)
    # dot on mirror line
    dot(ax, r_dot * np.cos(ang), r_dot * np.sin(ang))
ax.text(0, -1.15, 'hex., dots on\n3 mirror lines',
        ha='center', fontsize=7.5, style='italic')

# ---- D4: dots on all four mirror axes ----
ax = axes[8]
setup(ax, r'$D_4$')
sq_motif(ax)
r_mirror = 0.62
r_dot    = 0.40
mirror_angs = [0, np.pi/4, np.pi/2, 3*np.pi/4]
for ang in mirror_angs:
    ax.plot([ r_mirror*np.cos(ang),  -r_mirror*np.cos(ang)],
            [ r_mirror*np.sin(ang),  -r_mirror*np.sin(ang)],
            'r--', lw=0.9, alpha=0.6)
for ang in np.linspace(0, 2*np.pi, 9)[:-1]:   # 8 dots
    dot(ax, r_dot * np.cos(ang), r_dot * np.sin(ang), s=38)
ax.text(0, -1.15, 'sq., dots on\n4 mirror axes',
        ha='center', fontsize=7.5, style='italic')

# ---- D6: dots on all six mirror axes ----
ax = axes[9]
setup(ax, r'$D_6$')
hex_motif(ax)
r_mirror = 0.75
r_dot    = 0.50
for ang in np.linspace(0, np.pi, 6):   # 6 mirror lines
    ax.plot([-r_mirror*np.cos(ang),  r_mirror*np.cos(ang)],
            [-r_mirror*np.sin(ang),  r_mirror*np.sin(ang)],
            'r--', lw=0.9, alpha=0.55)
for ang in np.linspace(0, 2*np.pi, 7)[:-1]:   # 6 dots
    dot(ax, r_dot * np.cos(ang), r_dot * np.sin(ang), s=38)
ax.text(0, -1.15, 'hex., dots on\n6 mirror axes',
        ha='center', fontsize=7.5, style='italic')

# ------------------------------------------------------------------
# Layout and export
# ------------------------------------------------------------------

plt.tight_layout(pad=0.8, h_pad=1.8)

plt.savefig('Fig-1.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig('Fig-1.pdf',          bbox_inches='tight', facecolor='white')

print("Saved  Fig-1.png  (300 dpi)  and  Fig-1.pdf  (vector)")
plt.close()
