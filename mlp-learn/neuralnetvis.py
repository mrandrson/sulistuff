import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Load data
data = pd.read_csv("synthetic_data.csv")
x_data = data['x'].values
y_data = data['y'].values

## New ##
x_min, x_max = np.min(x_data), np.max(x_data)
y_min, y_max = np.min(y_data), np.max(y_data)

 
x_data_norm = (x_data - x_min) / (x_max - x_min)
x_fit = np.linspace(x_min, x_max, 300)
x_fit_norm = (x_fit - x_min) / (x_max - x_min)
##########

results = pd.read_csv("results.csv")
param_cols = [col for col in results.columns if col.startswith("p")]
params_over_time = results[param_cols].values
iterations = results['Iter'].values

def mlp_model(x, params, H1=30, H2=60, H3=30):
    x = np.array(x)
    offset = 0

    # Layer 1 (Input → H1)
    W1 = params[offset:offset+H1]
    b1 = params[offset+H1:offset+2*H1]
    h1 = np.tanh(np.outer(x, W1) + b1)
    offset += 2 * H1

    # Layer 2 (H1 → H2)
    W2 = params[offset:offset+H1*H2].reshape(H2, H1)
    b2 = params[offset+H1*H2:offset+H1*H2+H2]
    h2 = np.tanh(h1 @ W2.T + b2)
    offset += H1 * H2 + H2

    # Layer 3 (H2 → H3)
    W3 = params[offset:offset+H2*H3].reshape(H3, H2)
    b3 = params[offset+H2*H3:offset+H2*H3+H3]
    h3 = np.tanh(h2 @ W3.T + b3)
    offset += H2 * H3 + H3

    # Output layer (H3 → 1)
    W4 = params[offset:offset+H3]
    b4 = params[offset+H3]
    y = h3 @ W4 + b4

    return y


# Set up plot
fig, ax = plt.subplots(figsize=(8, 5))
x_fit = np.linspace(min(x_data), max(x_data), 300)
scatter = ax.scatter(x_data, y_data, color='gray', alpha=0.5, label="Noisy Data")
line, = ax.plot([], [], color='red', linewidth=2, label="MLP Fit")
title = ax.set_title("")
ax.set_xlim(min(x_data), max(x_data))
ax.set_ylim(min(y_data) - 0.5, max(y_data) + 0.5)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.legend()

# Interpolation setup
N_frames = 300
interp_params = np.linspace(0, len(params_over_time) - 1, N_frames)

def init():
    line.set_data([], [])
    title.set_text("")
    return line, title

def update(frame):
    lower_idx = int(np.floor(interp_params[frame]))
    upper_idx = min(lower_idx + 1, len(params_over_time) - 1)
    alpha = interp_params[frame] - lower_idx

    params = (1 - alpha) * params_over_time[lower_idx] + alpha * params_over_time[upper_idx]
    #y_fit = mlp_model(x_fit, params)
    ## New ##
    y_fit_norm = mlp_model(x_fit_norm, params)
    y_fit = y_fit_norm * (y_max - y_min) + y_min
    #########
    line.set_data(x_fit, y_fit)
    iter_number = int((1 - alpha) * iterations[lower_idx] + alpha * iterations[upper_idx])
    #title.set_text(f"Iteration {iter_number}")
    return line, title

ani = FuncAnimation(
    fig, update, frames=N_frames,
    init_func=init, blit=True, interval=75, repeat=False
)

#ani.save("/Users/richardanderson/adam_train.mp4")

plt.tight_layout()
plt.show()



