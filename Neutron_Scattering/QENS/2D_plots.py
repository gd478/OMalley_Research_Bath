import numpy as np
import matplotlib.pyplot as plt

# Define the full array of 18 Q values (in Å⁻¹)
all_Q_values = np.array([0.292721, 0.443938, 0.592860, 0.705534, 0.837721,
                         0.953758, 1.06934, 1.18212, 1.29340, 1.39237,
                         1.47924, 1.56702, 1.65104, 1.72433, 1.78835,
                         1.84431, 1.89158, 1.92991])
# Use only the first 6 Q values
Q_values = all_Q_values[:6]

# Load QENS data from file, skipping the header row.
data = np.loadtxt("qens_data.txt", skiprows=1)

# Extract energy transfer values (in MeV) from the first column.
x = data[:, 0]

# Number of Q channels to plot (first 6 Q channels)
num_channels = 6

# Create subplots: one row per Q channel.
fig, axs = plt.subplots(num_channels, 1, figsize=(8, 4 * num_channels), sharex=True)

# Loop over the first 6 Q channels
for i in range(num_channels):
    # Column indexing for each Q channel:
    # Measured data: column 1 + 6*i
    # Total fit: column 3 + 6*i
    # Lorentzian fit (f3 comp, column F): column 6 + 6*i
    y = data[:, 1 + i * 6]
    ftotal = data[:, 3 + i * 6]
    lorentzian = data[:, 6 + i * 6]
    
    ax = axs[i] if num_channels > 1 else axs
    # Plot measured data as small black crosses.
    ax.plot(x, y, 'kx', markersize=4, label="Data")
    # Plot the total fit as a thin black line.
    ax.plot(x, ftotal, 'k-', linewidth=1, label="Total Fit")
    # Plot the Lorentzian fit as a thin red line.
    ax.plot(x, lorentzian, 'r-', linewidth=1, label="Lorentzian Fit")
    
    ax.set_ylabel("Intensity (a.u.)")
    ax.legend(loc="best")

axs[-1].set_xlabel("Energy Transfer (MeV)")

plt.tight_layout()
plt.savefig("qens_plot.png", dpi=300)
plt.show()
