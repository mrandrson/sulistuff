from math import pi
import numpy as np
import matplotlib.pyplot as plt
import openmc
import openmc.model

fuel = openmc.Material(name='Fuel')
fuel.set_density('g/cm3', 10.5)
fuel.add_nuclide('U235', 4.6716e-02)
fuel.add_nuclide('U238', 2.8697e-01)
fuel.add_nuclide('O16',  5.0000e-01)
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

graphite = openmc.Material()
graphite.set_density('g/cm3', 1.1995)
graphite.add_element('C', 1.0)
graphite.add_s_alpha_beta('c_Graphite')

materials = openmc.Materials([fuel, buff, PyC1, PyC2, SiC, graphite])
materials.export_to_xml()

# --- Pebble generation ---
#pebble_centers = [[0, 0, 0], [6, 6, 6]]
pebble_centers = np.loadtxt('/Users/richardanderson/Downloads/pebble_centers_01.txt')[:, :3]
triso_template = np.loadtxt('/Users/richardanderson/Downloads/triso_centers.txt')[:, :3]
outer_radius = 4.25e-4  # cm
pebble_radius = 3.0  # cm
pebble_cells = []

for c in pebble_centers:
    # Create TRISO universe
    spheres = [openmc.Sphere(r=1e-4*r)
               for r in [215., 315., 350., 385.]]
    cells = [openmc.Cell(fill=fuel, region=-spheres[0]),
             openmc.Cell(fill=buff, region=+spheres[0] & -spheres[1]),
             openmc.Cell(fill=PyC1, region=+spheres[1] & -spheres[2]),
             openmc.Cell(fill=SiC, region=+spheres[2] & -spheres[3]),
             openmc.Cell(fill=PyC2, region=+spheres[3])]
    triso_univ = openmc.Universe(cells=cells)

    #region = -openmc.Sphere(r=3)
    region = -openmc.Sphere(x0=c[0], y0=c[1], z0=c[2], r=pebble_radius)
    outer_radius = 425.*1e-4
    # centers = openmc.model.pack_spheres(radius=outer_radius, region=region, pf=0.3, seed=124848351)
    centers = np.loadtxt('/Users/richardanderson/Downloads/triso_centers.txt')[:, :3]
    print(f"There are {len(centers)} trisos")
    shifted_centers = centers + np.array(c)
    trisos = [openmc.model.TRISO(outer_radius, triso_univ, center) for center in shifted_centers]
    #trisos = [openmc.model.TRISO(outer_radius, triso_univ, center) for center in centers]

    centers = np.vstack([triso.center for triso in trisos])

    box = openmc.Cell(region=region)
    lower_left, upper_right = box.region.bounding_box
    shape = (3, 3, 3)
    pitch = (upper_right - lower_left)/shape
    lattice = openmc.model.create_triso_lattice(
        trisos, lower_left, pitch, shape, graphite)

    box.fill = lattice
    pebble_cells.append(box)

universe = openmc.Universe(cells = pebble_cells)

geometry = openmc.Geometry(universe)
geometry.export_to_xml()

materials = list(geometry.get_all_materials().values())
openmc.Materials(materials).export_to_xml()

settings = openmc.Settings()
settings.run_mode = 'plot'
settings.export_to_xml()

# --- Plotting ---
shape = (120.0, 120.0, 120.0)
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

