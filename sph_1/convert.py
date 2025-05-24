import os
import meshio
import numpy as np

INPUT_FILE  = "post/sph_.dat"
OUTPUT_DIR  = "post_vtk"
os.makedirs(OUTPUT_DIR, exist_ok=True)

PEBBLE_RADIUS = 0.005  # [m]

'''
def write_vtk(timestep, atoms, cols):
    """Write one VTK file for this timestep."""
    # extract coordinates
    xi = cols.index('x')
    yi = cols.index('y')
    zi = cols.index('z')
    pts = atoms[:, [xi, yi, zi]]

    # extract types (cast to int)
    types = atoms[:, cols.index('type')].astype(int)

    # uniform radius
    radii = np.full(len(pts), PEBBLE_RADIUS)

    # build VTK mesh
    n_pts = len(pts)
    cells = [("vertex", np.arange(n_pts).reshape(-1, 1))]
    point_data = {"type": types, "radius": radii}

    mesh = meshio.Mesh(points=pts, cells=cells, point_data=point_data)
    outpath = os.path.join(OUTPUT_DIR, f"sph_{timestep}.vtk")
    meshio.write(outpath, mesh, file_format="vtk")
    print(f"Wrote {outpath}")
'''

def write_vtk(timestep, atoms, cols):
    xi, yi, zi = cols.index('x'), cols.index('y'), cols.index('z')
    pts   = atoms[:, [xi, yi, zi]]

    types = atoms[:, cols.index('type')].astype(int)
    radii = np.full(len(pts), PEBBLE_RADIUS)

    fxi, fyi, fzi = cols.index('fx'), cols.index('fy'), cols.index('fz')
    forces = atoms[:, [fxi, fyi, fzi]]

    cells = [("vertex", np.arange(len(pts)).reshape(-1,1))]
    point_data = {
      "type":   types,
      "radius": radii,
      "force" : forces
    }

    mesh = meshio.Mesh(points=pts, cells=cells, point_data=point_data)
    outpath = os.path.join(OUTPUT_DIR, f"sph_{timestep}.vtk")
    meshio.write(outpath, mesh, file_format="vtk")
    print(f"Wrote {outpath}")


cols      = None
atoms_buf = []
current_ts = None

with open(INPUT_FILE) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue

        # TIMESTEP header
        if line.startswith("ITEM: TIMESTEP"):
            if current_ts is not None and atoms_buf:
                atoms = np.array(atoms_buf, dtype=float)
                write_vtk(current_ts, atoms, cols)
                atoms_buf.clear()
            current_ts = f.readline().strip()
            continue

        # ATOMS header (with your long list of fields)
        if line.startswith("ITEM: ATOMS"):
            # split off the first two tokens and use the rest as column names
            cols = line.split()[2:]
            continue

        # data lines: only if we know the columns
        if cols is not None:
            tok = line.split()
            if len(tok) == len(cols):
                atoms_buf.append([float(x) for x in tok])

# flush last frame
if current_ts is not None and atoms_buf:
    atoms = np.array(atoms_buf, dtype=float)
    write_vtk(current_ts, atoms, cols)

