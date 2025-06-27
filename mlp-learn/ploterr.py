import pandas as pd
import matplotlib.pyplot as plt

# Load optimization results
results = pd.read_csv("results.csv")

# Extract iteration and error columns
iterations = results['Iter']
errors = results['Error']

# Plot
plt.figure(figsize=(8, 5))
plt.plot(iterations, errors, color='blue', linewidth=2)
plt.xlabel("Iteration")
plt.ylabel("Error (MSE)")
plt.yscale('log')
plt.xscale('log')
plt.title("MLP Training Error vs. Iteration")
plt.grid(True)
plt.tight_layout()
plt.show()

