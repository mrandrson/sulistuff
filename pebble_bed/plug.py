import numpy as np
import trimesh
import os

mesh_dir = os.path.join(os.path.dirname(__file__), "meshes")
stl_in   = os.path.join(mesh_dir, "GPBR_smooth_ascii_hollow.stl")

hollow_shell = trimesh.load(stl_in, process=False)

min_xy = hollow_shell.bounds[0, :2]
max_xy = hollow_shell.bounds[1, :2]
center_xy = (min_xy + max_xy) / 2

outlet_z = float(hollow_shell.bounds[0, 2])

d_inlet = 0.520
thickness = 0.05
r_outlet_inner = (d_inlet / 2) - thickness

plug_thickness = thickness

sections = 64

verts, faces = [], []

thetas = np.linspace(0, 2 * np.pi, sections, endpoint=False)
for ang in thetas:
    x = r_outlet_inner * np.cos(ang) + center_xy[0]
    y = r_outlet_inner * np.sin(ang) + center_xy[1]
    verts.append([x, y, outlet_z])
    verts.append([x, y, outlet_z - plug_thickness])

for i in range(sections):
    i0, i1 = 2 * i,     2 * i + 1
    j0, j1 = 2 * ((i + 1) % sections), 2 * ((i + 1) % sections) + 1
    faces.append([i0, j0, j1])
    faces.append([i0, j1, i1])

top_center = len(verts)
verts.append([center_xy[0], center_xy[1], outlet_z])
for i in range(sections):
    faces.append([top_center, 2*i, 2*((i+1)%sections)])

bottom_center = len(verts)
verts.append([center_xy[0], center_xy[1], outlet_z - plug_thickness])
for i in range(sections):
    faces.append([bottom_center, 2*((i+1)%sections) + 1, 2*i + 1])

plug = trimesh.Trimesh(vertices=np.array(verts), faces=np.array(faces), process=False)
stl_out = os.path.join(mesh_dir, "plug.stl")
plug.export(stl_out, file_type='stl')
print(f"Wrote properly centered, correctly sized plug: {stl_out}")

