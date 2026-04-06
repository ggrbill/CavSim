import sys
import csv
import numpy as np
import matplotlib.pyplot as plt

def plot_fields(csv_file, nv):
    """
    Reads a CSV file containing pressure field data and plots the pressure field.

    Args:
        csv_file (str): Path to the CSV file containing the data.
        nv (int): Number of divisions in the grid (assumes a square grid).

    Raises:
        ValueError: If the number of data points is not a perfect square.

    The CSV file is expected to have the following columns:
        - Column 1: x-coordinates
        - Column 2: y-coordinates
        - Column 3: Velocity in the x-direction (not used)
        - Column 4: Velocity in the y-direction (not used)
        - Column 5: Pressure values

    The function reshapes the data into a 2D grid and plots the pressure field
    using a color map. It also overlays the grid lines for better visualization.
    """
    # Read the CSV file
    x_vals = []
    y_vals = []
    p_vals = []

    with open(csv_file, 'r') as f:
        reader = csv.reader(f)
        next(reader)  # Skip the header
        for row in reader:
            x_vals.append(float(row[0]))
            y_vals.append(float(row[1]))
            p_vals.append(float(row[4]))  # Pressure is in the 5th column (index 4)

    num_data = len(x_vals)
    print(f"Data read: {num_data} points")
    print(f"Expected: {nv * nv} points (nv={nv})")

    if num_data != nv * nv:
        print(f"Warning: Number of points ({num_data}) does not match nv²={nv*nv}")
        # Attempt to determine the actual grid size
        import math
        real_nv = int(math.sqrt(num_data))
        if real_nv * real_nv == num_data:
            print(f"Correcting nv to {real_nv}")
            nv = real_nv
        else:
            raise ValueError(f"Number of points {num_data} is not a perfect square")

    # Reshape into 2D arrays (nv x nv)
    x = np.array(x_vals).reshape((nv, nv))
    y = np.array(y_vals).reshape((nv, nv))
    p = np.array(p_vals).reshape((nv, nv))

    # Calculate dx and dy (assuming uniform grid)
    if nv > 1:
        dx = x[0, 1] - x[0, 0]
        dy = y[1, 0] - y[0, 0]
    else:
        dx = 1.0  # Default value if nv=1
        dy = 1.0

    # Plot the pressure field
    fig, ax = plt.subplots(figsize=(8, 6))
    c = ax.pcolormesh(x, y, p, shading='auto', cmap='viridis')
    plt.colorbar(c, ax=ax, label='Pressure')

    # Add grid lines (black and thin)
    for i in range(nv + 1):
        ax.axvline(x=i * dx, color='black', linewidth=0.5)
    for j in range(nv + 1):
        ax.axhline(y=j * dy, color='black', linewidth=0.5)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('Pressure Field')
    ax.set_aspect('equal')  # Maintain square aspect ratio

    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Uso: python plot_pressure.py <arquivo_csv> <nv>")
        sys.exit(1)
    plot_fields(sys.argv[1], int(sys.argv[2]))