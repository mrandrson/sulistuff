'''
import numpy as np
import openmc

centers = np.loadtxt('/Users/richardanderson/Downloads/pebble_centers_02.txt')[:, :3]
radii = np.sqrt(centers[:, 0]**2 + centers[:, 1]**2)

print(f"There are {len(centers)}")
graphite = openmc.Material(name='Graphite')
graphite.add_element('C', 1.0)
graphite.set_density('g/cm3', 1.8)
materials = openmc.Materials([graphite])

sphere_cells = []
sphere_radius = 3.0
for i, (x, y, z) in enumerate(centers):
    sphere_surface = openmc.Sphere(x0=x, y0=y, z0=z, r=sphere_radius)
    sphere_cell = openmc.Cell(name=f'sphere_{i}', region=-sphere_surface)
    sphere_cell.fill = graphite
    sphere_cells.append(sphere_cell)

cyl = openmc.ZCylinder(r=50.0)
zmin = openmc.ZPlane(z0=-100.0)
zmax = openmc.ZPlane(z0=100.0)
cylinder_region = -cyl & +zmin & -zmax
outer_cell = openmc.Cell(name='outer_cylinder', region=cylinder_region)
outer_universe = openmc.Universe(cells=[outer_cell] + sphere_cells)
geometry = openmc.Geometry(root=outer_universe)
geometry.export_to_xml()
materials.export_to_xml()

settings = openmc.Settings()
settings.batches = 10
settings.inactive = 0
settings.particles = 1000
settings.run_mode = 'fixed source'
source = openmc.Source()
source.space = openmc.stats.Point((0, 0, 50))
settings.source = source
settings.export_to_xml()

shape = (100.0, 100.0, 100.0)
vox = openmc.Plot()
vox.type = 'voxel'
vox.filename = 'geom'
vox.width = shape
vox.pixels = (200, 200, 200)
vox.color_by = 'material'
plots = openmc.Plots([vox])
plots.export_to_xml()
openmc.plot_geometry(vox)
openmc.voxel_to_vtk('geom.h5', 'geom.vti')
'''

import numpy as np
import openmc
import time

# Load pebble centers
centers = np.loadtxt('/Users/richardanderson/Downloads/pebble_centers_02.txt')[:, :3]
sphere_radius = 3.0  # cm
print(f"There are {len(centers)} graphite pebbles.")

# --- Materials ---
fuel = openmc.Material(name='Fuel')
fuel.set_density('g/cm3', 10.5)
fuel.add_nuclide('U235', 4.6716e-02)
fuel.add_nuclide('U238', 2.8697e-01)
fuel.add_nuclide('O16', 5.0000e-01)
fuel.add_element('C', 1.6667e-01)

buff = openmc.Material(name='Buffer')
buff.set_density('g/cm3', 1.0)
buff.add_element('C', 1.0)
buff.add_s_alpha_beta('c_Graphite')

PyC1 = openmc.Material(name='PyC1')
PyC1.set_density('g/cm3', 1.9)
PyC1.add_element('C', 1.0)
PyC1.add_s_alpha_beta('c_Graphite')

PyC2 = openmc.Material(name='PyC2')
PyC2.set_density('g/cm3', 1.87)
PyC2.add_element('C', 1.0)
PyC2.add_s_alpha_beta('c_Graphite')

SiC = openmc.Material(name='SiC')
SiC.set_density('g/cm3', 3.2)
SiC.add_element('C', 0.5)
SiC.add_element('Si', 0.5)

graphite = openmc.Material(name='Graphite')
graphite.set_density('g/cm3', 1.8)
graphite.add_element('C', 1.0)
graphite.add_s_alpha_beta('c_Graphite')

materials = openmc.Materials([fuel, buff, PyC1, PyC2, SiC, graphite])
materials.export_to_xml()

# --- TRISO universe ---
spheres = [openmc.Sphere(r=1e-4 * r) for r in [215., 315., 350., 385.]]
cells = [
    openmc.Cell(fill=fuel, region=-spheres[0]),
    openmc.Cell(fill=buff, region=+spheres[0] & -spheres[1]),
    openmc.Cell(fill=PyC1, region=+spheres[1] & -spheres[2]),
    openmc.Cell(fill=SiC, region=+spheres[2] & -spheres[3]),
    openmc.Cell(fill=PyC2, region=+spheres[3]),
]
triso_univ = openmc.Universe(cells=cells)

# --- Load relative TRISO centers ---
rel_triso_centers = np.loadtxt('/Users/richardanderson/Downloads/triso_centers.txt')[:, :3]
# rel_triso_centers = rel_triso_centers[:2] ##Take only first two
outer_radius = 4.25e-4  # cm
triso_cells = []
start = time.time()

sphere_cells = []
for i, (x, y, z) in enumerate(centers):
    # Bounding region for pebble
    pebble_region = -openmc.Sphere(x0=x, y0=y, z0=z, r=sphere_radius)
    box_cell = openmc.Cell(region=pebble_region)

    # Define bounding box and shape of lattice
    lower_left = np.array([x - sphere_radius, y - sphere_radius, z - sphere_radius])
    upper_right = np.array([x + sphere_radius, y + sphere_radius, z + sphere_radius])
    shape = (int(2 * sphere_radius / 0.5),) * 3 
    pitch = (upper_right - lower_left) / shape

    # Place TRISOs at correct absolute positions
    triso_centers = rel_triso_centers + np.array([x, y, z])
    trisos = [openmc.model.TRISO(outer_radius, triso_univ, c) for c in triso_centers]
    print(trisos)
    # Create lattice and assign to pebble
    lattice = openmc.model.create_triso_lattice(trisos, lower_left, pitch, shape, graphite)
    box_cell.fill = lattice
    sphere_cells.append(box_cell)

end = time.time()
print(f"TRISOs packed into all pebbles in {end - start:.2f} seconds.")

# --- Outer cylinder (container)
cyl = openmc.ZCylinder(r=50.0)
zmin = openmc.ZPlane(z0=-100.0)
zmax = openmc.ZPlane(z0=100.0)
cylinder_region = -cyl & +zmin & -zmax
outer_cell = openmc.Cell(name='outer_cylinder', region=cylinder_region)

# --- Combine geometry ---
outer_universe = openmc.Universe(cells=[outer_cell] + sphere_cells)
geometry = openmc.Geometry(root=outer_universe)
geometry.export_to_xml()

# --- Settings ---
settings = openmc.Settings()
settings.batches = 10
settings.inactive = 0
settings.particles = 1000
settings.run_mode = 'fixed source'
source = openmc.Source()
source.space = openmc.stats.Point((0, 0, 50))
settings.source = source
settings.export_to_xml()

'''
shape = (100.0, 100.0, 100.0)
vox = openmc.Plot()
vox.type = 'voxel'
vox.filename = 'geom'
vox.width = shape
vox.pixels = (200, 200, 200)
vox.color_by = 'material'
plots = openmc.Plots([vox])
plots.export_to_xml()
openmc.plot_geometry(vox)
openmc.voxel_to_vtk('geom.h5', 'geom.vti')
'''
