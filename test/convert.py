import os
import glob
import numpy as np
import meshio

input_dir = "post"
output_dir = "post_vtk"

os.makedirs(output_dir, exist_ok=True)

for filepath in sorted(glob.glob(os.path.join(input_dir, "dump*.newModels"))):
    basename = os.path.splitext(os.path.basename(filepath))[0]
    output_path = os.path.join(output_dir, f"{basename}.vtk")
    
    with open(filepath, 'r') as f:
        lines = f.readlines()

    start_idx = None
    for idx, line in enumerate(lines):
        if line.startswith("ITEM: ATOMS"):
            cols = line.strip().split()[2:]
            start_idx = idx + 1
            break
    if start_idx is None:
        print(f"Skipping {filepath}: no 'ITEM: ATOMS' header found.")
        continue

    data_lines = lines[start_idx:]
    if not data_lines or all(not ln.strip() for ln in data_lines):
        print(f"Skipping {filepath}: no data lines found.")
        continue

    try:
        data = np.loadtxt(data_lines)
    except Exception as e:
        print(f"Error reading numeric data in {filepath}: {e}")
        continue

    if data.ndim == 1:
        print(f"Skipping {filepath}: only one data row found, expected multiple.")
        continue

    col_indices = {col: i for i, col in enumerate(cols)}

    points = data[:, [col_indices['x'], col_indices['y'], col_indices['z']]]
    radii = data[:, col_indices['radius']]
    types = data[:, col_indices['type']].astype(int)

    cells = [("vertex", np.arange(len(points)).reshape(-1, 1))]
    point_data = {
        "radius": radii,
        "type": types
    }

    mesh = meshio.Mesh(points=points, cells=cells, point_data=point_data)
    try:
        meshio.write(output_path, mesh, file_format="vtk")
        print(f"Converted {filepath} → {output_path}")
    except Exception as e:
        print(f"Failed to write {output_path}: {e}")

