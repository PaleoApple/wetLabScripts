import numpy as np
import matplotlib.pyplot as plt

# -------------------------------
# Load CSV
# -------------------------------
rt = []
intensity_matrix = []

with open("SPPS_spectaP.csv") as f:
    lines = f.readlines()

# First row = m/z bins
first_line = lines[0].strip().rstrip(">").split(",")
mz_bins = np.array([float(x) for x in first_line])

# Remaining rows
for line in lines[1:]:
    line = line.strip().rstrip(">")
    fields = line.split(",")
    if len(fields) < 2:
        continue
    try:
        rt_val = float(fields[0])
        intensities = []
        for x in fields[1:]:
            try:
                intensities.append(float(x))
            except ValueError:
                intensities.append(0.0)
        if len(intensities) < len(mz_bins):
            intensities += [0.0]*(len(mz_bins) - len(intensities))
        elif len(intensities) > len(mz_bins):
            intensities = intensities[:len(mz_bins)]
        rt.append(rt_val)
        intensity_matrix.append(intensities)
    except ValueError:
        continue

rt = np.array(rt)
intensity_matrix = np.array(intensity_matrix)

# -------------------------------
# Mask RT and m/z ranges first
# -------------------------------
rt_min, rt_max = 240.0, 1400.0  # seconds
mask = (rt >= rt_min) & (rt <= rt_max)
rt = rt[mask]
intensity_matrix = intensity_matrix[mask, :]

mz_min, mz_max = 600, 1200
mz_mask = (mz_bins >= mz_min) & (mz_bins <= mz_max)
mz_bins = mz_bins[mz_mask]
intensity_matrix = intensity_matrix[:, mz_mask]

print("RT shape:", rt.shape)
print("Intensity matrix shape:", intensity_matrix.shape)
print("mz_bins shape:", mz_bins.shape)

# -------------------------------
# TIC Plot
# -------------------------------
tic = intensity_matrix.sum(axis=1)
'''
plt.plot(rt/60, tic)
plt.xlabel("Retention Time (min)")
plt.ylabel("Total Ion Intensity")
plt.title("TIC")
plt.show()
'''
# -------------------------------
# Load peptide library
# -------------------------------
library = {}
with open("permutations_output.dat") as f:
    next(f)  # skip header
    for line in f:
        name, mass = line.split()
        library[name] = float(mass)

PROTON = 1.007276
def mz_theoretical(M, z):
    return (M + z*PROTON)/z

tolerance = 0.8  # Da
top_N = 2         # top N peaks per scan
min_score = 0.6   # Δm/z scoring threshold
min_consecutive = 1  # minimum scans for a valid peptide

# -------------------------------
# Scan-based peak assignment
# -------------------------------
results = []

for i, rt_val in enumerate(rt):
    spectrum = intensity_matrix[i, :]
    
    # scan-specific noise
    scan_median = np.median(spectrum)
    scan_std = np.std(spectrum)
    noise_mask = spectrum > (scan_median + 2*scan_std)
    
    if not np.any(noise_mask):
        results.append((rt_val, "?", None, None, None))
        continue
    
    valid_idx = np.where(noise_mask)[0]
    # top-N peaks among significant ones
    top_idx = valid_idx[np.argsort(spectrum[valid_idx])[-top_N:][::-1]]
    top_mz = mz_bins[top_idx]
    top_int = spectrum[top_idx]

    best_peptide = None
    best_z = None
    best_mz = None
    best_int = 0

    for mz_obs, intensity in zip(top_mz, top_int):
        for peptide, mass in library.items():
            for z in [2,1]:
                target_mz = mz_theoretical(mass, z)
                score = 1 - abs(mz_obs - target_mz)/tolerance
                if score >= min_score and intensity > best_int:
                    best_peptide = peptide
                    best_z = z
                    best_mz = mz_obs
                    best_int = intensity

    if best_peptide:
        results.append((rt_val/60, best_peptide, best_z, best_mz, best_int))
    else:
        results.append((rt_val/60, "?", None, None, None))

# -------------------------------
# Merge consecutive same-peptide assignments
# -------------------------------
merged_results = []
i = 0
while i < len(results):
    current_peptide = results[i][1]
    if current_peptide == "?":
        i += 1
        continue
    # count consecutive appearances
    count = 1
    for j in range(i+1, len(results)):
        if results[j][1] == current_peptide:
            count += 1
        else:
            break
    if count >= min_consecutive:
        group = results[i:i+count]
        apex = max(group, key=lambda x: x[4] if x[4] is not None else 0)
        merged_results.append(apex)
    i += count

# -------------------------------
# Print final peptide assignments
# -------------------------------
print("RT (s) | Peptide | z | Observed m/z | Intensity")
for r in merged_results:
    print(r)

plt.figure(figsize=(12,5))
plt.plot(rt/60, tic, color='blue', label='TIC')
plt.xlabel("Retention Time (min)")
plt.ylabel("Total Ion Intensity")
plt.title("TIC with Peptide Apex Annotations")

# annotate merged peptide apexes using TIC intensity
for r in merged_results:
    rt_apex = r[0]*60  # apex in seconds
    peptide_name = r[1]

    # find closest RT index in rt array
    idx = np.argmin(np.abs(rt - rt_apex))
    tic_value = tic[idx]

    # plot on TIC
    plt.scatter(rt[idx]/60, tic_value, color='red', s=50)
    plt.text(rt[idx]/60, tic_value*1.05, peptide_name, rotation=45,
             fontsize=12, ha='left', va='bottom', color='black')

plt.tight_layout()
plt.legend()
plt.show()

