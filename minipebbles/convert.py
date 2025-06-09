'''
#!/usr/bin/env python3
import os
import glob
import numpy as np
import meshio

input_dir = "post"
output_dir = "post_vtk"
os.makedirs(output_dir, exist_ok=True)

for src in sorted(glob.glob(os.path.join(input_dir, "packing_*.dat"))):
    basename = os.path.splitext(os.path.basename(src))[0]
    dst = os.path.join(output_dir, f"{basename}.vtk")

    # Read all lines
    with open(src) as f:
        lines = f.readlines()

    # Find the line index just after "ITEM: ATOMS"
    start_idx = next(
        (i+1 for i, line in enumerate(lines) if line.strip().startswith("ITEM: ATOMS")),
        None
    )
    if start_idx is None:
        print(f"Skipping {src}: no ITEM: ATOMS header found")
        continue

    # Load the numeric block
    data = np.loadtxt(lines[start_idx:])

    # Columns: id=0, type=1, x=2, y=3, z=4, …, radius=last
    pts = data[:, 2:5]
    radii = data[:, -1]
    types = data[:, 1].astype(int)

    # Build a VTK point-cloud (vertices)
    cells = [("vertex", np.arange(len(pts)).reshape(-1, 1))]
    point_data = {
        "radius": radii,
        "type": types
    }

    mesh = meshio.Mesh(
        points=pts,
        cells=cells,
        point_data=point_data
    )

    meshio.write(dst, mesh, file_format="vtk")
    print(f"Converted {src} → {dst}")
'''

import os
import meshio

INPUT_FILE  = "post/chute.dat"
OUTPUT_DIR  = "post_vtk"

os.makedirs(OUTPUT_DIR, exist_ok=True)

def write_vtk(timestep, atoms, cols):
    """Write one VTK file for this timestep."""
    pts   = atoms[:, cols.index('x'):cols.index('x')+3]
    types = atoms[:, cols.index('type')].astype(int)
    radii = atoms[:, cols.index('radius')]

    n_pts = len(pts)
    cells = [("vertex", np.arange(n_pts).reshape(-1, 1))]
    point_data = {"type": types, "radius": radii}

    mesh = meshio.Mesh(points=pts, cells=cells, point_data=point_data)
    outpath = os.path.join(OUTPUT_DIR, f"chute_particles_{timestep}.vtk")
    meshio.write(outpath, mesh, file_format="vtk")
    print(f"Wrote {outpath}")

import numpy as np

cols      = None
atoms_buf = []
current_ts = None

with open(INPUT_FILE) as f:
    for line in f:
        tok = line.strip().split()
        if not tok:
            continue

        if tok[0] == "ITEM:" and tok[1] == "TIMESTEP":
            if current_ts is not None and atoms_buf:
                atoms = np.array(atoms_buf, dtype=float)
                write_vtk(current_ts, atoms, cols)
                atoms_buf.clear()

            current_ts = f.readline().strip()

        elif tok[0] == "ITEM:" and tok[1] == "ATOMS":
            cols = tok[2:]

        elif cols is not None and len(tok) == len(cols):
            atoms_buf.append([float(x) for x in tok])

if current_ts is not None and atoms_buf:
    atoms = np.array(atoms_buf, dtype=float)
    write_vtk(current_ts, atoms, cols)

