import openmc
import numpy as np
import matplotlib.pyplot as plt

uo2 = openmc.Material(name='UO2 fuel')
uo2.add_element('U', 1.0, enrichment=5.0)
uo2.add_element('O', 2.0)
uo2.set_density('g/cm3', 10.5)

graphite = openmc.Material(name='Graphite')
graphite.add_element('C', 1.0)
graphite.set_density('g/cm3', 1.7)

materials = openmc.Materials([uo2, graphite])

fuel_sphere = openmc.Sphere(r=2.0, boundary_type='vacuum')

fuel_region = -fuel_sphere

fuel_cell = openmc.Cell(name='Fuel pebble', region=fuel_region, fill=uo2)
root_universe = openmc.Universe(cells=[fuel_cell])

geometry = openmc.Geometry(root_universe)
materials.export_to_xml()
geometry.export_to_xml()

source = openmc.Source()
source.space = openmc.stats.Point((0, 0, 0))
source.angle = openmc.stats.Isotropic()
source.energy = openmc.stats.Watt(a=0.988e6, b=2.249e-6)

settings = openmc.Settings()
settings.batches = 100
settings.inactive = 10
settings.particles = 10000
settings.source = source
settings.output = {'tallies': True}
settings.export_to_xml()

energy_bins = openmc.mgxs.GROUP_STRUCTURES['CASMO-70']
energy_filter = openmc.EnergyFilter(energy_bins)

neutron_flux = openmc.Tally(name='Neutron flux')
neutron_flux.filters = [energy_filter]
neutron_flux.scores = ['flux']

gamma_spectrum = openmc.Tally(name='Gamma spectrum')
gamma_spectrum.filters = [energy_filter]
gamma_spectrum.scores = ['flux']

tallies = openmc.Tallies([neutron_flux, gamma_spectrum])
tallies.export_to_xml()

openmc.run()

sp = openmc.StatePoint('statepoint.100.h5')

tally = sp.get_tally(name='Neutron flux')
neutron_flux_data = tally.get_pandas_dataframe()

tally_gamma = sp.get_tally(name='Gamma spectrum')
gamma_flux_data = tally_gamma.get_pandas_dataframe()

plt.figure(figsize=(8,6))
plt.loglog(neutron_flux_data['energy low [eV]'], neutron_flux_data['mean'], label='Neutron Flux')
plt.xlabel('Energy (eV)')
plt.ylabel('Flux')
plt.legend()
plt.show()

plt.figure(figsize=(8,6))
plt.loglog(gamma_flux_data['energy low [eV]'], gamma_flux_data['mean'], label='Gamma Spectrum', color='red')
plt.xlabel('Energy (eV)')
plt.ylabel('Intensity')
plt.legend()
plt.show()

