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

