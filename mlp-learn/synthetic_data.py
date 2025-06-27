import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# 1. Generate synthetic data
np.random.seed(0)
x = np.linspace(0.1, 10, 100)  # Avoid x=0 for log safety
y = (
    1.2*np.exp(x-x**2) * np.sin(2 * np.pi * x / 5) +                          # Periodic feature
    np.exp(-(x - 3.0)**2 / (2 * 0.2**2)) * 2.0 +               # Sharp peak
    np.exp(-(x - 7.0)**2 / (2 * 0.5**2)) * 1.2 +               # Wider peak
    np.random.normal(0, 0.2, size=x.shape)                    # Stochastic noise
)

# 2. Plot for inspection
plt.plot(x, y, label='y = true function with noise')
plt.xlabel("x")
plt.ylabel("y")
plt.title("Synthetic Data: Peaks + Periodic + Noise")
plt.legend()
plt.grid()
plt.show()

# 3. Save to CSV for C++ use
df = pd.DataFrame({
    "x": x,
    "y": y
})
df.to_csv("synthetic_data.csv", index=False)

print("Saved synthetic_data.csv with", len(x), "points.")


