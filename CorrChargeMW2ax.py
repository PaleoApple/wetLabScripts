import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ----------------------------
# INPUT: Sequence
# ----------------------------
J = "XYZ"
sequence = J  # your residue string
positions = np.arange(1, len(sequence) + 1)  # biological numbering

# ----------------------------
# PROPERTY MAPS
# ----------------------------
charge_map = {
    'K': 1.0, 'R': 1.0, 'H': 0.04,
    'D': -1.0, 'E': -1.0,
    'C': -0.1, 'Y': 0.0
}

mass_map = {
    'G': 57.051, 'A': 71.078, 'S': 87.077, 'P': 97.115,
    'V': 99.131, 'T': 101.104, 'C': 103.143,
    'I': 113.158, 'L': 113.158,
    'N': 114.103, 'D': 115.087,
    'Q': 128.129, 'K': 128.172,
    'E': 129.114, 'M': 131.196,
    'H': 137.139, 'F': 147.174,
    'R': 156.186, 'Y': 163.173, 'W': 186.210
}

# ----------------------------
# FUNCTIONS
# ----------------------------
def compute_property_trace(sequence, property_map):
    return np.array([property_map.get(res, 0.0) for res in sequence])

def rolling_trace(values, window=5):
    return pd.Series(abs(values)).rolling(window=window, center=True).mean().values

def zscore_trace(values):
    return (values - np.nanmean(values)) / np.nanstd(values)

def residue_color(aa, i, gray_after=40):
    if i > gray_after:
        return 'gray'
    # coloring by chemistry
    if aa in ['K','R']:
        return 'blue'
    elif aa in ['D','E']:
        return 'red'
    elif aa == 'H':
        return 'purple'
    else:
        return 'black'

def normalize_minus1_to1(arr):
    arr_min = np.min(arr)
    arr_max = np.max(arr)
    # Avoid division by zero if all values are equal
    if arr_max - arr_min == 0:
        return np.zeros_like(arr)
    return (arr - arr_min) / (arr_max - arr_min)

# Apply to your rolling traces
# ----------------------------
# CALCULATE TRACES
# ----------------------------
charge_vals = compute_property_trace(sequence, charge_map)
mass_vals   = compute_property_trace(sequence, mass_map)

window = 5
half_window = window // 2

# Compute rolling traces WITHOUT NaNs at the edges
charge_roll_full = rolling_trace(charge_vals, window=window)
mass_roll_full   = rolling_trace(mass_vals, window=window)

# Trim edges
charge_roll = charge_roll_full[half_window:-half_window]
mass_roll   = mass_roll_full[half_window:-half_window]

# Corresponding positions for plotting
positions_roll = positions[half_window:-half_window]

# Normalize
charge_roll = normalize_minus1_to1(charge_roll)
mass_rolln   = normalize_minus1_to1(mass_roll)

print(charge_roll)

#charge_z = zscore_trace(charge_roll)
#mass_z   = zscore_trace(mass_roll)

delta = np.abs(charge_roll - mass_rolln)
threshold = 0.2
mask = delta > threshold

# Create masked arrays to keep only regions above threshold
charge_masked = np.ma.masked_where(delta <= threshold, charge_roll)
mass_masked   = np.ma.masked_where(delta <= threshold, mass_rolln)

# ----------------------------
# PLOT
# ----------------------------
# ----------------------------
# PLOT
# ----------------------------
fig, ax1 = plt.subplots(figsize=(9,3))

# LEFT AXIS — Normalized Charge
ax1.plot(positions_roll, charge_roll,
         color='black', linewidth=2, label='Abs. Charge')
ax1.set_xlabel("Residue", labelpad=15)
ax1.set_ylabel(f"Normalized 5-residue\nAbs. Charge & Mass", color='black')
ax1.tick_params(axis='y', labelcolor='black')

# RIGHT AXIS — Raw Rolling Mass
ax2 = ax1.twinx()

ax2.plot(positions_roll, mass_roll,
         color='orange', linewidth=2, label='Mass (Da/Norm)')
ax2.set_ylabel("Rolling Mass Average (Da)", color='orange')
ax2.tick_params(axis='y', labelcolor='orange')

# ----------------------------
# Difference shading
# (must normalize mass temporarily for comparison only)
# ----------------------------

mass_norm_for_delta = normalize_minus1_to1(mass_roll)
delta = np.abs(charge_roll - mass_norm_for_delta)
threshold = 0.2

charge_masked = np.ma.masked_where(delta <= threshold, charge_roll)
mass_masked   = np.ma.masked_where(delta <= threshold, mass_rolln)

ax1.fill_between(
    positions_roll,
    charge_masked,
    mass_masked,
    color='purple',
    alpha=0.3,
    interpolate=True,
    label='Norm Diff > 0.2'
)

# ----------------------------
# Residue annotations
# ----------------------------
ymin, ymax = ax1.get_ylim()
y_text = ymin - 0.15*(ymax - ymin)

for i, aa in enumerate(sequence, start=1):
    ax1.text(
        i,
        y_text,
        aa,
        ha='center',
        va='top',
        fontsize=8,
        family='monospace',
        color=residue_color(aa, i)
)

# Combine legends
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='lower left')

plt.xlim(1, min(positions[-1], 43))
plt.subplots_adjust(bottom=0.25)
plt.suptitle("Overlay of Absolute Charge & Raw Mass (norm) Rolling Average for VL42C",
             y=0.94)

plt.tight_layout()
plt.savefig('OverlayChargeMW_VL42C', dpi=600)
plt.show()
