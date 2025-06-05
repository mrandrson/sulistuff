import numpy as np
import trimesh

factor    = 0.1

d_inlet   = factor * 0.520
h_inlet   = factor * 1.851
h_taper   = factor * 0.55438
d_main    = factor * 2.400
h_main    = factor * 9.388
h_outlet  = factor * 0.200
thickness = factor * 0.05

def tube_segment(r0, r1, height, z_offset, sections=64):
    theta = np.linspace(0, 2*np.pi, sections, endpoint=False)
    verts, faces = [], []
    for ang in theta:
        verts.append([r0*np.cos(ang), r0*np.sin(ang), 0])
        verts.append([r1*np.cos(ang), r1*np.sin(ang), height])
    for i in range(sections):
        i0, i1 = 2*i, 2*i+1
        j0, j1 = 2*((i+1)%sections), 2*((i+1)%sections)+1
        faces.append([i0, j0, j1])
        faces.append([i0, j1, i1])
    mesh = trimesh.Trimesh(np.array(verts), np.array(faces), process=False)
    mesh.apply_translation([0, 0, z_offset])
    return mesh

r_in_o = d_inlet/2
r_main_o = d_main/2
r_in_i = r_in_o - thickness
r_main_i = r_main_o - thickness

outer_inlet = tube_segment(r_in_o, r_in_o, h_inlet, 0)
outer_taper = tube_segment(r_in_o, r_main_o, h_taper, h_inlet)
outer_main  = tube_segment(r_main_o, r_main_o, h_main, h_inlet+h_taper)
outer_out   = tube_segment(r_main_o, r_main_o, h_outlet, h_inlet+h_taper+h_main)

inner_inlet = tube_segment(r_in_i, r_in_i, h_inlet, 0)
inner_taper = tube_segment(r_in_i, r_main_i, h_taper, h_inlet)
inner_main  = tube_segment(r_main_i, r_main_i, h_main, h_inlet+h_taper)
inner_out   = tube_segment(r_main_i, r_main_i, h_outlet, h_inlet+h_taper+h_main)

# reverse inner faces
for seg in [inner_inlet, inner_taper, inner_main, inner_out]:
    seg.faces = seg.faces[:, ::-1]

# assemble hollow shell
hollow_shell = trimesh.util.concatenate([
    outer_inlet, outer_taper, outer_main, outer_out,
    inner_inlet, inner_taper, inner_main, inner_out
])

# ----------------------------------------------------------
# Positioning block: move mouth to z=0.15 and center to (-0.2,0)
# ----------------------------------------------------------
# compute current bounds
bounds = hollow_shell.bounds  # shape (2,3): [[min x,y,z], [max x,y,z]]
current_top_z    = bounds[1, 2]
current_center_xy = (bounds[0, :2] + bounds[1, :2]) / 2

# desired placement
desired_top_z    = 0.15
desired_center_xy = np.array([0, 0])  # = [-0.2, 0.0]

# compute translation vector
dx, dy = desired_center_xy - current_center_xy
dz     = desired_top_z    - current_top_z

# apply
hollow_shell.apply_translation([dx, dy, dz])

# export
output_path ='/Users/richardanderson/workdir/suliwork/minipebbles/meshes/GPBR_smooth_ascii_hollow.stl'
hollow_shell.export(output_path, file_type='stl')
print(f"ASCII STL written to {output_path}")

