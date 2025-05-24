#!/usr/bin/env python
import openmc
import openmc.lib
import matplotlib.pyplot as plt

material_list = []

uo2 = openmc.Material(material_id=1, name='UO2 fuel at 2.4% wt enrichment')
uo2.set_density('g/cm3', 10.29769)
uo2.add_element('U', 1., enrichment=2.4)
uo2.add_element('O', 2.)
material_list.append(uo2)

helium = openmc.Material(material_id=2, name='Helium for gap')
helium.set_density('g/cm3', 0.001598)
helium.add_element('He', 2.4044e-4)
material_list.append(helium)

zircaloy = openmc.Material(material_id=3, name='Zircaloy 4')
zircaloy.set_density('g/cm3', 6.55)
zircaloy.add_element('Sn', 0.014, 'wo')
zircaloy.add_element('Fe', 0.00165, 'wo')
zircaloy.add_element('Cr', 0.001, 'wo')
zircaloy.add_element('Zr', 0.98335, 'wo')
material_list.append(zircaloy)

for i in range(4, 14):
    water = openmc.Material(material_id=i)
    water.set_density('g/cm3', 0.7)
    water.add_element('H', 2.0)
    water.add_element('O', 1.0)
    water.add_s_alpha_beta('c_H_in_H2O')
    material_list.append(water)
    
materials_file = openmc.Materials(material_list)
materials_file.export_to_xml()

pitch = 1.25984
fuel_or = openmc.ZCylinder(r=0.39218)
clad_ir = openmc.ZCylinder(r=0.40005)
clad_or = openmc.ZCylinder(r=0.4572)
left = openmc.XPlane(x0=-pitch/2)
right = openmc.XPlane(x0=pitch/2)
back = openmc.YPlane(y0=-pitch/2)
front = openmc.YPlane(y0=pitch/2)
z = [0., 30., 60., 90., 120., 150., 180., 210., 240., 270., 300.]
z_list = [openmc.ZPlane(z0=z_i) for z_i in z]

left.boundary_type = 'reflective'
right.boundary_type = 'reflective'
front.boundary_type = 'reflective'
back.boundary_type = 'reflective'
z_list[0].boundary_type = 'vacuum'
z_list[-1].boundary_type = 'vacuum'

fuel_list = []
gap_list = []
clad_list = []
water_list = []
for i in range(1, 11):
    fuel_list.append(openmc.Cell(cell_id=i))
    gap_list.append(openmc.Cell(cell_id=i+10))
    clad_list.append(openmc.Cell(cell_id=i+20))
    water_list.append(openmc.Cell(cell_id=i+30))
    
for j, fuels in enumerate(fuel_list):
    fuels.region = -fuel_or & +z_list[j] & -z_list[j+1]
    fuels.fill = uo2
    fuels.temperature = 800.

for j, gaps in enumerate(gap_list):
    gaps.region = +fuel_or & -clad_ir & +z_list[j] & -z_list[j+1]
    gaps.fill = helium
    gaps.temperature = 700.

for j, clads in enumerate(clad_list):
    clads.region = +clad_ir & -clad_or & +z_list[j] & -z_list[j+1]
    clads.fill = zircaloy
    clads.temperature = 600.

for j, waters in enumerate(water_list):
    waters.region = +clad_or & +left & -right & +back & -front & +z_list[j] & -z_list[j+1]
    waters.fill = material_list[j+3]
    waters.temperature = 500.

root = openmc.Universe(name='root universe')
root.add_cells(fuel_list)
root.add_cells(gap_list)
root.add_cells(clad_list)
root.add_cells(water_list)
geometry_file = openmc.Geometry(root)
geometry_file.export_to_xml()

cell_filter = openmc.CellFilter(fuel_list)
t = openmc.Tally(tally_id=1)
t.filters.append(cell_filter)
t.scores = ['fission-q-recoverable']
tallies = openmc.Tallies([t])
tallies.export_to_xml()

root.plot(basis='yz', width=[2, 10], color_by='material', origin=[0., 0., 150.], pixels=[400, 400])

#plt.show()

lower_left = [-0.62992, -pitch/2, 0]
upper_right = [+0.62992, +pitch/2, +300]
uniform_dist = openmc.stats.Box(lower_left, upper_right, only_fissionable=True)

settings_file = openmc.Settings()
settings_file.batches = 100
settings_file.inactive = 10
settings_file.particles = 10000
settings_file.temperature = {'multipole': True, 'method': 'interpolation', 'range': [290, 2500]}
settings_file.source = openmc.source.IndependentSource(space=uniform_dist)
settings_file.export_to_xml()

openmc.lib.init()
openmc.lib.simulation_init()

for _ in range(14):
    openmc.lib.next_batch()

t = openmc.lib.tallies[1]
print(t.mean)

print("fuel temperature is: ")
print(openmc.lib.cells[5].get_temperature())
print("gap temperature is: ")
print(openmc.lib.cells[15].get_temperature())
print("clad temperature is: ")
print(openmc.lib.cells[25].get_temperature())
print("water temperature is: ")
print(openmc.lib.cells[35].get_temperature())

for i in range(1, 11):
    temp = 900.0
    openmc.lib.cells[i].set_temperature(temp)

print("fuel temperature is: ")
print(openmc.lib.cells[5].get_temperature())

for i in range(4, 14):
    density = 0.65
    openmc.lib.materials[i].set_density(density, units='g/cm3')

for _ in range(14):
    openmc.lib.next_batch()

openmc.lib.simulation_finalize()
openmc.lib.finalize()

import openmc
import openmc.lib
import matplotlib.pyplot as plt
import numpy as np

axial_positions = np.linspace(0, 300, 10)  # Midpoints of axial zones

# Ensure simulation is initialized
if not openmc.lib.is_initialized:
    openmc.lib.init()
    openmc.lib.simulation_init()

# Extract cell IDs from OpenMC's in-memory model
cell_ids = [cell.id for cell in openmc.lib.cells.values()]

# Create tally with correct filter
tally = openmc.lib.Tally(24)
cell_filter = openmc.lib.TallyFilter(openmc.lib.TallyFilter.CELL, cell_ids)  # Corrected line
tally.filters = [cell_filter]
tally.scores = ['flux']
tally.activate()  # Register the tally

# Run simulation to collect tally results
for _ in range(14):
    openmc.lib.next_batch()

# Retrieve tally results
flux_tally = openmc.lib.tallies[24]
flux_means = flux_tally.mean

# Plot neutron flux distribution
plt.figure(figsize=(8,6))
plt.plot(axial_positions, flux_means, marker='o', linestyle='-', color='b', label="Neutron Flux")
plt.xlabel("Axial Position (cm)")
plt.ylabel("Flux")
plt.title("Axial Neutron Flux Distribution")
plt.legend()
plt.grid()
plt.show()

openmc.lib.simulation_finalize()
openmc.lib.finalize()

initial_temperatures = [800.0] * 10
updated_temperatures = [openmc.lib.cells[i].get_temperature() for i in range(1, 11)]

plt.figure(figsize=(8,6))
plt.plot(axial_positions, initial_temperatures, marker='o', linestyle='--', color='r', label="Initial Temperature")
plt.plot(axial_positions, updated_temperatures, marker='o', linestyle='-', color='b', label="Updated Temperature")
plt.xlabel("Axial Position (cm)")
plt.ylabel("Temperature (K)")
plt.title("Axial Fuel Temperature Changes")
plt.legend()
plt.grid()
plt.show()

power_deposition = flux_means * 200  # Approximate assumption

plt.figure(figsize=(8,6))
plt.plot(axial_positions, power_deposition, marker='o', linestyle='-', color='g', label="Fission Power")
plt.xlabel("Axial Position (cm)")
plt.ylabel("Power Deposition (arbitrary units)")
plt.title("Axial Fission Power Deposition")
plt.legend()
plt.grid()
plt.show()

y_positions = np.linspace(-pitch/2, pitch/2, 10)
flux_grid = np.random.rand(len(y_positions), len(axial_positions))  # Placeholder for real flux data

plt.figure(figsize=(8,6))
plt.contourf(y_positions, axial_positions, flux_grid.T, levels=20, cmap='viridis')
plt.colorbar(label="Neutron Flux")
plt.xlabel("Radial Position (cm)")
plt.ylabel("Axial Position (cm)")
plt.title("2D Neutron Flux Contour")
plt.show()

initial_densities = [0.7] * 10
updated_densities = [openmc.lib.materials[i].get_density() for i in range(4, 14)]

plt.figure(figsize=(8,6))
plt.plot(axial_positions, initial_densities, marker='o', linestyle='--', color='r', label="Initial Density")
plt.plot(axial_positions, updated_densities, marker='o', linestyle='-', color='b', label="Updated Density")
plt.xlabel("Axial Position (cm)")
plt.ylabel("Density (g/cm³)")
plt.title("Axial Water Density Changes")
plt.legend()
plt.grid()
plt.show()

