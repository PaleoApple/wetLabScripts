import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# List of files corresponding to the 4 conditions
files = [
    "DHP 0mM NaCl.csv",
    "DHP 50mM NaCl.csv",
    "DHP 150mM NaCl.csv",
    "DHP 300mM NaCl.csv"
]

ymin = -0.05
ymax = 1

def process_file_data(file_path):
    if not os.path.exists(file_path):
        print(f"Warning: File not found -> {file_path}")
        return None, None
        
    skip_rows = 0
    with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
        for i, line in enumerate(f):
            if "Temperature (°C),Abs" in line:
                skip_rows = i
                break

    df = pd.read_csv(file_path, skiprows=skip_rows)
    if df.columns[0] == 'Unnamed: 0' or df.columns[0] == '':
        df = df.iloc[:, 1:]

    sample_names = []
    with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if "SAMPLES" in line and (i + 1) < len(lines):
                raw_cols = [c.strip() for c in lines[i + 1].split(',')]
                sample_names = [col for col in raw_cols if col not in ['Name', ''] and 'Unnamed' not in col]
                break
                
    return df, sample_names


# Initialize the 2x2 grid layout
fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharex=True, sharey=True)
axes = axes.flatten() 

for idx, file_path in enumerate(files):
    ax = axes[idx]
    df, sample_names = process_file_data(file_path)
    
    if df is None:
        ax.text(0.5, 0.5, f"Missing Data:\n{file_path}", ha='center', va='center', fontsize=12, color='gray')
        ax.set_title(file_path.replace(".csv", ""), fontsize=12, fontweight='bold')
        continue

    num_pairs = df.shape[1] // 2

    for i in range(0, num_pairs):
        temp_col_idx = i * 2
        abs_col_idx = i * 2 + 1
        
        temp_vals = pd.to_numeric(df.iloc[:, temp_col_idx], errors='coerce')
        abs_vals = pd.to_numeric(df.iloc[:, abs_col_idx], errors='coerce')
        
        plot_data = pd.concat([temp_vals, abs_vals], axis=1).dropna()
        if plot_data.empty:
            continue
            
        # Ignore row 0 initialization spikes
        x = plot_data.iloc[1:, 0].values
        y = plot_data.iloc[1:, 1].values
        
        base_label = sample_names[i] if i < len(sample_names) else f"Sample {i+1}"
        
        # --- DETECT STAGE FLIP AND STRIP COOLING ---
        diffs = np.diff(x)
        flip_indices = np.where(np.diff(np.sign(diffs)))[0]
        
        if len(flip_indices) > 0:
            flip_point = flip_indices[0] + 1
            
            # Stage 1: Heating only (Cooling block completely stripped out)
            x_heat, y_heat = x[:flip_point], y[:flip_point]
            y_heat_zeroed = y_heat - y_heat[0] 
            # FIXED: Label clean with no "(H)" suffix
            ax.plot(x_heat, y_heat_zeroed, label=base_label, linewidth=1.2)
                
        else:
            y_zeroed = y - y[0]
            ax.plot(x, y_zeroed, label=base_label, linewidth=1.2)

    condition_title = file_path.replace(".csv", "")
    ax.set_title(condition_title, fontsize=12, fontweight='bold')
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.set_ylim(ymin, ymax)
    
    ax.legend(loc='upper left', fontsize=8, framealpha=0.6, ncol=2)

fig.suptitle(r"Temperature vs. Baseline-Corrected Absorbance ($\Delta$Abs) for DHP across NaCl Conditions", fontsize=16, fontweight='bold', y=0.98)

for ax in axes[2:]: 
    ax.set_xlabel("Temperature (°C)", fontsize=12)
for ax in [axes[0], axes[2]]: 
    ax.set_ylabel(r"$\Delta$ Absorbance (AU)", fontsize=12)

plt.tight_layout()
plt.savefig('DHPmultiRunUvVisPlot.png', dpi=400)
plt.show()
