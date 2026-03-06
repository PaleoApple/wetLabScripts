import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import matplotlib as mpl
import numpy as np

# -------------------------------
# Publication-quality formatting
# -------------------------------
mpl.rcParams.update({
    "savefig.dpi": 600,
    "font.size": 10,
    "axes.labelsize": 10,
    "axes.titlesize": 12,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9
})

# ===============================
# ========= FILE PATHS ==========
# ===============================

samples = [
    {
        "name": "hLAD C1 WT (Isabelle)",
        "color": "tab:blue",
        "hplc": "C1_TIC.csv",
        "mz": "C1_mz.csv",
        "deconv": "C1_dec.csv"
    },
    {
        "name": "hLAD C1 R219M (Allison)",
        "color": "tab:red",
        "hplc": "C1R219M_TIC.csv",
        "mz": "C1R219_mz.csv",
        "deconv": "C1R219M_dec.csv"
    },
    {
        "name": "hLAD C2 WT (Oscar)",
        "color": "tab:green",
        "hplc": "ovrerl/cleavedC2tic_OW.csv",
        "mz": "ovrerl/cleavedC2mz_OW.csv",
        "deconv": "ovrerl/cleavedC2decon_OW.csv"
    }
]

# ===============================
# ========= PARAMETERS ==========
# ===============================

detect_peaks = True

peak_params = {
    "hplc": {"height": 1e11, "distance": 10},
    "mz": {"height": 1e8, "distance": 5},
    "deconv": {"height": 4e8, "distance": 500}
}

normalize = False  # Max-normalize each trace if True
max_labels = 10    # Maximum number of peak annotations per sample

# ===============================
# ========= FUNCTIONS ===========
# ===============================


'''
def load_xy(file):
    df = pd.read_csv(file)
    x = df.iloc[:, 0]
    y = df.iloc[:, 1]
    return x, y
'''

def load_xy(file, rt_in_minutes=False):
    df = pd.read_csv(file)
    x = df.iloc[:, 0]
    y = df.iloc[:, 1]
    if rt_in_minutes:
        x = x / 60  # convert seconds to minutes
    return x, y

def maybe_normalize(y):
    return y / y.max() if normalize else y

def overlay_panel_multi(ax, samples, key, label, xlabel, ylabel, peak_key, rt_in_minutes=False):
    for s_idx, sample in enumerate(samples):
        x, y = load_xy(sample[key], rt_in_minutes=rt_in_minutes)
        y = maybe_normalize(y)

        ax.plot(
            x, y,
            label=sample["name"],
            linewidth=1.2,
            color=sample["color"]
        )

        if detect_peaks:
            params = peak_params[peak_key]
            peaks, _ = find_peaks(
                y,
                height=params["height"],
                distance=params["distance"]
            )

            # Peak markers
            ax.plot(
                x.iloc[peaks],
                y.iloc[peaks],
                "o",
                markersize=4,
                color=sample["color"]
            )

            # Limit annotations
            peaks_sorted = peaks[np.argsort(y.iloc[peaks])][::-1][:max_labels]

            # Axis-based offset
            x_min, x_max = ax.get_xlim()
            y_min, y_max = ax.get_ylim()
            vertical_base = (y_max - y_min) * 0.03
            horizontal_offset = (x_max - x_min) * 0.005

            # Stagger vertically by sample index
            for idx, i in enumerate(peaks_sorted):
                ax.text(
                    x.iloc[i] + horizontal_offset,
                    y.iloc[i] + vertical_base * (1 + idx % 2 + s_idx),
                    f"{x.iloc[i]:.2f}",
                    fontsize=9,
                    color=sample["color"],
                    ha="left",
                    bbox=dict(
                        boxstyle="round,pad=0.2",
                        facecolor="white",
                        edgecolor=sample["color"],
                        linewidth=0.6,
                        alpha=0.9
                    )
                )

    ax.set_title(label)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


# HPLC Overlay with RT in minutes
'''
def overlay_panel_rt_minutes(ax, file_A, file_B, label, xlabel, ylabel, peak_key):
    # Only convert RT for this plot
    xA, yA = load_xy(file_A, rt_in_minutes=True)
    xB, yB = load_xy(file_B, rt_in_minutes=True)

    yA = maybe_normalize(yA)
    yB = maybe_normalize(yB)

    color_A = "tab:blue"
    color_B = "tab:red"

    ax.plot(xA, yA, label="Sample A", linewidth=1.2, color=color_A)
    ax.plot(xB, yB, label="Sample B", linewidth=1.2, color=color_B)

    if detect_peaks:
        params = peak_params[peak_key]
        peaks_A, _ = find_peaks(yA, height=params["height"], distance=params["distance"])
        peaks_B, _ = find_peaks(yB, height=params["height"], distance=params["distance"])
        ax.plot(xA.iloc[peaks_A], yA.iloc[peaks_A], "o", markersize=4, color=color_A)
        ax.plot(xB.iloc[peaks_B], yB.iloc[peaks_B], "o", markersize=4, color=color_B)

    ax.set_title(label)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
'''

def overlay_panel_rt_minutes(ax, file_A, file_B, label, xlabel, ylabel, peak_key):
    # Convert RT to minutes
    xA, yA = load_xy(file_A, rt_in_minutes=True)
    xB, yB = load_xy(file_B, rt_in_minutes=True)

    yA = maybe_normalize(yA)
    yB = maybe_normalize(yB)

    color_A = "tab:blue"
    color_B = "tab:red"

    ax.plot(xA, yA, label="hLAD C1 WT (Isabelle)", linewidth=1.2, color=color_A)
    ax.plot(xB, yB, label="hLAD C1 R219M (Allison)", linewidth=1.2, color=color_B)

    if detect_peaks:
        params = peak_params[peak_key]
        peaks_A, _ = find_peaks(yA, height=params["height"], distance=params["distance"])
        peaks_B, _ = find_peaks(yB, height=params["height"], distance=params["distance"])

        # Plot markers
        ax.plot(xA.iloc[peaks_A], yA.iloc[peaks_A], "o", markersize=4, color=color_A)
        ax.plot(xB.iloc[peaks_B], yB.iloc[peaks_B], "o", markersize=4, color=color_B)

        # Axis range for offsets
        x_min, x_max = ax.get_xlim()
        y_min, y_max = ax.get_ylim()
        y_range = y_max - y_min
        vertical_base = y_range * 0.03
        horizontal_offset = (x_max - x_min) * 0.005  # small +x shift

        # Annotate Sample A
        for idx, i in enumerate(peaks_A[:10]):
            ax.text(
                xA.iloc[i] + horizontal_offset,
                yA.iloc[i] + vertical_base * (1 + idx % 2),
                f"{xA.iloc[i]:.2f}",
                fontsize=10,
                color=color_A,
                ha="left",
                bbox=dict(
                    boxstyle="round,pad=0.2",
                    facecolor="white",
                    edgecolor=color_A,
                    linewidth=0.6,
                    alpha=0.9
                )
            )

        # Annotate Sample B
        for idx, i in enumerate(peaks_B[:10]):
            ax.text(
                xB.iloc[i] + horizontal_offset,
                yB.iloc[i] + vertical_base * (2 + idx % 2),
                f"{xB.iloc[i]:.2f}",
                fontsize=10,
                color=color_B,
                ha="left",
                bbox=dict(
                    boxstyle="round,pad=0.2",
                    facecolor="white",
                    edgecolor=color_B,
                    linewidth=0.6,
                    alpha=0.9
                )
            )

    # Formatting
    ax.set_title(label)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)



# ===============================
# ========= MAIN FIGURE =========
# ===============================

fig, axes = plt.subplots(3, 1, figsize=(16, 10))

# HPLC
overlay_panel_multi(
    axes[0],
    samples,
    key="hplc",
    label="HPLC Overlay",
    xlabel="Retention Time (min)",
    ylabel="TIC",
    peak_key="hplc",
    rt_in_minutes=True
)

# m/z
overlay_panel_multi(
    axes[1],
    samples,
    key="mz",
    label="m/z Scan Overlay",
    xlabel="m/z",
    ylabel="Intensity",
    peak_key="mz"
)

# Deconvolved
overlay_panel_multi(
    axes[2],
    samples,
    key="deconv",
    label="Deconvolved Spectrum Overlay",
    xlabel="Mass (Da)",
    ylabel="Intensity",
    peak_key="deconv"
)

axes[2].set_ylim(0,4e9)

axes[0].legend()
plt.tight_layout()
plt.savefig("overlay_comparison_3samples.png", bbox_inches="tight")
plt.show()
# Optional save:

