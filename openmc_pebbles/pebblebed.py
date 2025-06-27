import math
from math import pi
import numpy as np
import openmc
import openmc.model
import openmc.deplete
import glob
import warnings
import os
import tracemalloc
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

## Define Materials ##
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
######################

packing = True 

if packing:
    threading.Thread(target=monitor, daemon=True).start()
    spheres = [openmc.Sphere(r=1e-4 * r) for r in [215., 315., 350., 385.]]
    cells = [
        openmc.Cell(fill=fuel, region=-spheres[0]),
        openmc.Cell(fill=buff, region=+spheres[0] & -spheres[1]),
        openmc.Cell(fill=PyC1, region=+spheres[1] & -spheres[2]),
        openmc.Cell(fill=SiC, region=+spheres[2] & -spheres[3]),
        openmc.Cell(fill=PyC2, region=+spheres[3]),
    ]
    triso_univ = openmc.Universe(cells=cells)

    radius = 10.0 #was 0.5
    height = 20.0 #was 1.0

    '''    
    outer_cyl = openmc.ZCylinder(r=radius, boundary_type='vacuum')
    bottom = openmc.ZPlane(z0=-height / 2, boundary_type='vacuum')
    top = openmc.ZPlane(z0=+height / 2, boundary_type='vacuum')
    region = -outer_cyl & +bottom & -top
    '''    

    region = -openmc.Sphere(r=3)

    outer_radius = 4.25e-4

    print("Packing TRISO particles.")
    start = time.time()
    #centers = openmc.model.pack_spheres(radius=outer_radius, region=region, pf=0.3, seed=124848351)
    center = np.loadtxt('/Users/richardanderson/Downloads/triso_centers.txt')
    centers = center[:, :3]
    print(centers)
    trisos = [openmc.model.TRISO(outer_radius, triso_univ, c) for c in centers]
    end = time.time()
    print(f"Packing complete. {len(centers)} TRISOs packed in {end - start:.2f} seconds.")
    box = openmc.Cell(region=region)
    lower_left, upper_right = box.region.bounding_box
    shape = (int(np.floor(2*radius/0.5)), int(np.floor(2*radius/0.5)), int(np.floor(height/0.5)))
    # ^ was (3, 3, 3), denominator greater than 0.4, less than 1. Smaller is more fine, slower, more is faster. 
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

    settings.particles = 50000
    settings.batches = 100
    settings.inactive = 20
    settings.verbosity = 10

    tallies = openmc.Tallies()
    mesh = openmc.RegularMesh()
    ll, ur = box.region.bounding_box
    mesh.lower_left  = ll
    mesh.upper_right = ur          
    mesh_pitch = 0.5
    bbox_size = upper_right - lower_left  
    mesh_dim = tuple(int(np.floor((bbox_size[i]) / mesh_pitch)) for i in range(3))
    mesh.dimension = mesh_dim
    mesh_filter = openmc.MeshFilter(mesh)

    tally = openmc.Tally(name='flux')
    tally.filters = [mesh_filter]
    tally.scores = ['flux', 'fission']
    tallies.append(tally)
    tallies.export_to_xml()

    model = openmc.Model(
        geometry,
        openmc.Materials(list(geometry.get_all_materials().values())),
        settings,
    )

    r_kernel = 0.025
    n_triso = 12000
    V_kernel = (4.0 / 3.0) * math.pi * r_kernel**3
    fuel.volume = len(trisos)*n_triso * V_kernel
    total_volume = len(trisos) * n_triso * V_kernel
    print(f"Fuel volume: {total_volume:.3e} cm^3")
    print(f"Number of TRISOs packed: {len(trisos)}")

    model.export_to_xml()

    '''     
    vox = openmc.Plot()
    vox.type       = 'voxel'
    vox.filename   = 'geom'
    vox.width = shape
    vox.pixels     = (200,200,200)
    vox.color_by   = 'material'
    plots = openmc.Plots([vox])
    plots.export_to_xml()
    openmc.plot_geometry(vox)
    openmc.voxel_to_vtk('geom.h5', 'geom.vti')
    '''    

    #openmc.run(tracks=True)
else:
    print("Skip Packing.")

plot = openmc.Plot()
plot.basis = 'xz'
plot.origin = (0.0, 0.0, 0.0)
plot.width = (20., 20.)
plot.pixels = (100, 100)
plots = openmc.Plots([plot])
plots.export_to_xml()
openmc.plot_geometry()


'''
height = 105.0
radius = 55.0
shape = (int(np.floor(2*radius/0.5)), int(np.floor(2*radius/0.5)), int(np.floor(height/0.5)))
vox = openmc.Plot()
vox.type       = 'voxel'
vox.filename   = 'geom'
vox.width = shape
vox.pixels     = (50,50,50)
vox.color_by   = 'material'
plots = openmc.Plots([vox])
plots.export_to_xml()
openmc.plot_geometry(vox)
openmc.voxel_to_vtk('geom.h5', 'geom.vti')
'''
#openmc.run()
#print(tracemalloc.get_traced_memory())

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

    # model.export_to_xml()

    model = openmc.model.Model.from_xml()
    chain = "/Users/richardanderson/openmc_data/chain_jeff33.xml"
    op = openmc.deplete.CoupledOperator(model, chain_file = chain)

    total_days = 100
    n_steps    = 5

    dt = total_days / n_steps

    timesteps = [dt] * n_steps
    power = 250
    openmc.deplete.CECMIntegrator(op, timesteps, power, timestep_units='d').integrate()
'''
