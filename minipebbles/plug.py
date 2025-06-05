import numpy as np
import trimesh
import os

# assume this script lives one level above "meshes/" folder
mesh_dir = os.path.join(os.path.dirname(__file__), "meshes")
stl_in   = os.path.join(mesh_dir, "GPBR_smooth_ascii_hollow.stl")

# load your translated chute
hollow_shell = trimesh.load(stl_in, process=False)

# compute mesh center in XY
min_xy = hollow_shell.bounds[0, :2]
max_xy = hollow_shell.bounds[1, :2]
center_xy = (min_xy + max_xy) / 2

# outlet is min z of the shell
outlet_z = float(hollow_shell.bounds[0, 2])

# inner radius at inlet (use d_inlet/2 - thickness)
d_inlet = 0.520
thickness = 0.05
r_outlet_inner = (d_inlet / 2) - thickness

# plug thickness (same as chute wall)
plug_thickness = thickness

# segmentation
sections = 64

verts, faces = [], []

# build side surface (with XY center offset)
thetas = np.linspace(0, 2 * np.pi, sections, endpoint=False)
for ang in thetas:
    x = r_outlet_inner * np.cos(ang) + center_xy[0]
    y = r_outlet_inner * np.sin(ang) + center_xy[1]
    verts.append([x, y, outlet_z])
    verts.append([x, y, outlet_z - plug_thickness])

# side faces
for i in range(sections):
    i0, i1 = 2 * i,     2 * i + 1
    j0, j1 = 2 * ((i + 1) % sections), 2 * ((i + 1) % sections) + 1
    faces.append([i0, j0, j1])
    faces.append([i0, j1, i1])

# top cap
top_center = len(verts)
verts.append([center_xy[0], center_xy[1], outlet_z])
for i in range(sections):
    faces.append([top_center, 2*i, 2*((i+1)%sections)])

# bottom cap
bottom_center = len(verts)
verts.append([center_xy[0], center_xy[1], outlet_z - plug_thickness])
for i in range(sections):
    faces.append([bottom_center, 2*((i+1)%sections) + 1, 2*i + 1])

# export sealed cylinder plug
plug = trimesh.Trimesh(vertices=np.array(verts), faces=np.array(faces), process=False)
stl_out = os.path.join(mesh_dir, "plug.stl")
plug.export(stl_out, file_type='stl')
print(f"Wrote properly centered, correctly sized plug: {stl_out}")

