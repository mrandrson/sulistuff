import math
from math import pi
import numpy as np
import openmc
import openmc.model
import openmc.deplete
import glob
import warnings
import os

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

graphite = openmc.Material()
graphite.set_density('g/cm3', 1.1995)
graphite.add_element('C', 1.0)
graphite.add_s_alpha_beta('c_Graphite')


spheres = [openmc.Sphere(r=1e-4 * r) for r in [215., 315., 350., 385.]]
cells = [
    openmc.Cell(fill=fuel, region=-spheres[0]),
    openmc.Cell(fill=buff, region=+spheres[0] & -spheres[1]),
    openmc.Cell(fill=PyC1, region=+spheres[1] & -spheres[2]),
    openmc.Cell(fill=SiC, region=+spheres[2] & -spheres[3]),
    openmc.Cell(fill=PyC2, region=+spheres[3]),
]
triso_univ = openmc.Universe(cells=cells)

radius = 0.5
height = 1.0
outer_cyl = openmc.ZCylinder(r=radius, boundary_type='vacuum')
bottom = openmc.ZPlane(z0=-height / 2, boundary_type='vacuum')
top = openmc.ZPlane(z0=+height / 2, boundary_type='vacuum')
region = -outer_cyl & +bottom & -top

outer_radius = 425.0e-4
centers = openmc.model.pack_spheres(radius=outer_radius, region=region, pf=0.3, seed=124848351)
trisos = [openmc.model.TRISO(outer_radius, triso_univ, c) for c in centers]

box = openmc.Cell(region=region)
lower_left, upper_right = box.region.bounding_box
shape = (3, 3, 3)
pitch = (upper_right - lower_left) / shape
lattice = openmc.model.create_triso_lattice(trisos, lower_left, pitch, shape, graphite)
box.fill = lattice

universe = openmc.Universe(cells=[box])
geometry = openmc.Geometry(universe)
fuel.depletable = True

materials = list(geometry.get_all_materials().values())
openmc.Materials(materials).export_to_xml()

settings = openmc.Settings()
settings.run_mode = 'eigenvalue'
settings.batches = 100
settings.inactive = 10
settings.particles = int(5e3)
settings.track = [(1, 2, 4)]
settings.output = {'tallies': False}

model = openmc.Model(
    geometry,
    openmc.Materials(list(geometry.get_all_materials().values())),
    settings,
)

r_kernel = 0.025
n_triso = 12000
V_kernel = (4.0 / 3.0) * math.pi * r_kernel**3
fuel.volume = len(trisos)*n_triso * V_kernel
model.export_to_xml()

vox = openmc.Plot()
vox.type       = 'voxel'
vox.filename   = 'geom'
vox.width      = (1.0, 1.0, 1.0)
vox.pixels     = (400, 400, 400)
vox.color_by   = 'material'
plots = openmc.Plots([vox])
plots.export_to_xml()
openmc.plot_geometry(vox)

openmc.voxel_to_vtk('geom.h5', 'geom.vti')

if __name__ == "__main__":
    for pat in ["openmc_simulation_*.h5",
                "statepoint.*.h5",
                "summary.h5",
                "depletion_results.h5",
                "tallies.out"]:
        for f in glob.glob(pat):
            try:
                os.remove(f)
            except OSError:
                pass

    model.export_to_xml()
    chain = "/Users/richardanderson/openmc_data/chain_jeff33.xml"
    openmc.run(tracks=True)
    '''
    op = openmc.deplete.CoupledOperator(model, chain_file = chain)

    total_days = 100
    n_steps    = 5

    dt = total_days / n_steps

    timesteps = [dt] * n_steps
    power = 250
    openmc.deplete.CECMIntegrator(op, timesteps, power, timestep_units='d').integrate()
    '''
