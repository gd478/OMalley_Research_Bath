import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.ticker as ticker

# Define the full array of 18 Q values (in Å⁻¹)
all_Q_values = np.array([0.292721, 0.443938, 0.592860, 0.705534, 0.837721,
                         0.953758, 1.06934, 1.18212, 1.29340, 1.39237,
                         1.47924, 1.56702, 1.65104, 1.72433, 1.78835,
                         1.84431, 1.89158, 1.92991])
# Use only the first 7 Q values
Q_values = all_Q_values[:7]

# Load QENS data from file, skipping the header row.
data = np.loadtxt("qens_data.txt", skiprows=1)

# Extract energy transfer values (in MeV) from the first column.
energy = data[:, 0]

# Number of Q channels to plot (first 7 channels)
num_channels = 7

# Create the 3D plot
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

# Loop over the first 7 Q channels
for i in range(num_channels):
    # For each Q channel, assume the following column structure:
    # - Measured data: column 1 + 6*i
    # - Total fit: column 3 + 6*i (plotted in red)
    # - Lorentzian fit: column 6 + 6*i (plotted in blue)
    y_measured = data[:, 1 + i * 6]
    total_fit  = data[:, 3 + i * 6]
    lorentzian = data[:, 6 + i * 6]
    
    # Create a constant Q array for the current channel
    current_Q = Q_values[i]
    Q_const = np.full_like(energy, current_Q)
    
    # Plot measured data as small black crosses.
    ax.scatter(energy, Q_const, y_measured, color='k', marker='x', s=20,
               label="Data" if i == 0 else "")
    # Plot total fit as a thin red line.
    ax.plot(energy, Q_const, total_fit, 'r-', linewidth=1,
            label="Total Fit" if i == 0 else "")
    # Plot Lorentzian fit as a thin blue line.
    ax.plot(energy, Q_const, lorentzian, 'b-', linewidth=1,
            label="Lorentzian Fit" if i == 0 else "")

# Set axis labels.
ax.set_xlabel("Energy Transfer (MeV)")
ax.set_ylabel("Q (Å⁻¹)")
ax.set_zlabel("Intensity (a.u.)")

# Reverse the Q axis (highest Q at the bottom)
ax.set_ylim(Q_values.max(), Q_values.min())

# Customize energy axis ticks: increments of 0.02, formatted to 1 significant figure.
ax.xaxis.set_major_locator(ticker.MultipleLocator(0.02))
ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: format(x, '.1g')))

# Set the background panes to light grey
ax.xaxis.pane.set_facecolor('lightgrey')
ax.yaxis.pane.set_facecolor('lightgrey')
ax.zaxis.pane.set_facecolor('lightgrey')

# Optionally, remove grid lines if not desired:
ax.grid(False)

# Set the plot box aspect ratio to be equal (cubed shape)
ax.set_box_aspect((1,1,1))

ax.legend(loc="best")
plt.tight_layout()
plt.savefig("qens_3d_plot.png", dpi=300)
plt.show()

