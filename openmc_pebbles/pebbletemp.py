import math
from math import pi
import numpy as np
import openmc
import openmc.model
import openmc.deplete
import glob
import warnings
import os
import psutil
import time
import threading

def monitor():
    pid = os.getpid()
    proc = psutil.Process(pid)
    with open("cpu_usage_log.txt", "w") as f:
        f.write("Time(s), CPU%, Memory(MB)\n")
        while True:
            cpu = proc.cpu_percent(interval=1)
            mem = proc.memory_info().rss / (1024 * 1024)
            timestamp = time.time()
            f.write(f"{timestamp}, {cpu}, {mem}\n")
            f.flush()

packing = True 

if packing:

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

    def pack_pebble(triso_univ, pebble_radius=3.0, triso_radius=0.0425, pf=0.3, seed=42):
        print("Packing Pebbles")
        pebble_sphere = openmc.Sphere(r=pebble_radius, boundary_type='transmission')
        region = -pebble_sphere

        centers = openmc.model.pack_spheres(radius=triso_radius, region=region, pf=pf, seed=seed)
        trisos = [openmc.model.TRISO(triso_radius, triso_univ, c) for c in centers]
        print(f"There are {len(trisos)} trisos.")

        pebble_cell = openmc.Cell(region=region)
        lower_left, upper_right = pebble_cell.region.bounding_box
        print(f"Lower Left : {lower_left}, Upper Right : {upper_right}")
        shape = (3, 3, 3)
        pitch = (upper_right - lower_left) / shape
        lattice = openmc.model.create_triso_lattice(trisos, lower_left, pitch, shape, graphite)
        pebble_cell.fill = lattice
        print(f"Lattice : {lattice}")

        pebble_universe = openmc.Universe(cells=[pebble_cell])
        
        return pebble_universe, len(trisos)

    def pack_core(triso_univ, core_radius, height, pf=0.1, seed=123):
        print("Packing Core")
        outer_cyl = openmc.ZCylinder(r=core_radius, boundary_type='vacuum')
        bottom = openmc.ZPlane(z0=-height/2, boundary_type='vacuum')
        top = openmc.ZPlane(z0=height/2, boundary_type='vacuum')
        region = -outer_cyl & +bottom & -top

        pebble_radius = 3.0
        centers = openmc.model.pack_spheres(
            radius=pebble_radius, region=region, pf=pf, seed=seed
        )

        filled_pebbles = []
        for i, c in enumerate(centers):
            pebble, _ = pack_pebble(triso_univ, seed=seed + i) 

            translated_pebble = openmc.Cell(
                region=+openmc.Sphere(x0=c[0], y0=c[1], z0=c[2], r=pebble_radius),
                fill=pebble
            )
            filled_pebbles.append(translated_pebble)
            print(f"There are {len(filled_pebbles)} pebbles.")
        return openmc.Universe(cells=filled_pebbles)

    spheres = [openmc.Sphere(r=1e-4 * r) for r in [215., 315., 350., 385.]]
    cells = [
        openmc.Cell(fill=fuel, region=-spheres[0]),
        openmc.Cell(fill=buff, region=+spheres[0] & -spheres[1]),
        openmc.Cell(fill=PyC1, region=+spheres[1] & -spheres[2]),
        openmc.Cell(fill=SiC, region=+spheres[2] & -spheres[3]),
        openmc.Cell(fill=PyC2, region=+spheres[3]),
    ]
    triso_univ = openmc.Universe(cells=cells)

    core_universe = pack_core(triso_univ, 10.0, 20.0, pf=0.3, seed=999)

    geometry = openmc.Geometry(core_universe)
    fuel.depletable = True

    materials = list(geometry.get_all_materials().values())
    openmc.Materials(materials).export_to_xml()
    geometry.export_to_xml()

    settings = openmc.Settings()
    settings.run_mode = 'eigenvalue'
    settings.batches = 100
    settings.inactive = 10
    settings.particles = int(1e4)
    settings.output = {'tallies': False}
    settings.export_to_xml()

    r_kernel = 0.025
    V_kernel = (4.0 / 3.0) * math.pi * r_kernel**3

    n_pebbles = len(core_universe.get_all_cells())
    _, n_trisos_per_pebble = pack_pebble(triso_univ, seed=42)
    fuel.volume = n_pebbles * n_trisos_per_pebble * V_kernel

    print(f"Total TRISOs estimated: {n_pebbles * n_trisos_per_pebble}")
    print(f"Fuel kernel volume: {fuel.volume:.3f} cm³")

    model = openmc.Model(geometry, openmc.Materials(materials), settings)
    model.export_to_xml()

model = openmc.model.Model.from_xml()

'''
vox = openmc.Plot()
vox.type = 'voxel'
vox.filename = 'geom'
vox.width = (30.0, 30.0, 30.0)
vox.pixels = (100, 100, 100)
vox.color_by = 'material'

plots = openmc.Plots([vox])
plots.export_to_xml()

openmc.plot_geometry()
openmc.voxel_to_vtk('geom.h5', 'geom.vti')
print("Exported geom.vti")
'''
'''
plot = openmc.Plot()
plot.basis = 'yz'
plot.origin = (0, 0, 0)
plot.width = (10.0, 10.0)
plot.pixels = (1200, 1200)
plot.color_by = 'material'

plots = openmc.Plots([plot])
plots.export_to_xml()
openmc.plot_geometry()
'''

'''
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

    model = openmc.model.Model.from_xml()
    chain = "/Users/richardanderson/openmc_data/chain_jeff33.xml"
    threading.Thread(target=monitor, daemon=True).start()
    #openmc.run(tracks=True)
    
    op = openmc.deplete.CoupledOperator(model, chain_file = chain)

    total_days = 100
    n_steps    = 5

    dt = total_days / n_steps

    timesteps = [dt] * n_steps
    power = 250
    openmc.deplete.CECMIntegrator(op, timesteps, power, timestep_units='d').integrate()
    
'''
