import openmc
import openmc.deplete
import numpy as np
import os

def run_burnup_simulation(enrichment, power, burn_times, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    uo2 = openmc.Material(name='UO2 fuel')
    uo2.add_element('U', 1.0, enrichment=enrichment)
    uo2.add_element('O', 2.0)
    uo2.set_density('g/cm3', 10.5)
    uo2.volume = 33.51

    graphite = openmc.Material(name='Graphite')
    graphite.add_element('C', 1.0)
    graphite.set_density('g/cm3', 1.7)

    materials = openmc.Materials([uo2, graphite])

    fuel_sphere = openmc.Sphere(r=2.0, boundary_type='vacuum')
    fuel_cell = openmc.Cell(name='Fuel pebble', region=-fuel_sphere, fill=uo2)
    root_universe = openmc.Universe(cells=[fuel_cell])

    geometry = openmc.Geometry(root_universe)

    source = openmc.IndependentSource()
    source.space = openmc.stats.Point((0, 0, 0))
    source.angle = openmc.stats.Isotropic()
    source.energy = openmc.stats.Watt(a=0.988e6, b=2.249e-6)

    settings = openmc.Settings()
    settings.batches = 100
    settings.inactive = 10
    settings.particles = 10000
    settings.source = source

    energy_bins = openmc.mgxs.GROUP_STRUCTURES['CASMO-70']
    energy_filter = openmc.EnergyFilter(energy_bins)

    neutron_flux = openmc.Tally(name='Neutron flux')
    neutron_flux.filters = [energy_filter]
    neutron_flux.scores = ['flux']

    gamma_spectrum = openmc.Tally(name='Gamma spectrum')
    gamma_spectrum.filters = [energy_filter]
    gamma_spectrum.scores = ['flux']

    tallies = openmc.Tallies([neutron_flux, gamma_spectrum])

    model = openmc.Model(geometry=geometry, materials=materials, settings=settings, tallies=tallies)

    #chain_file = "chain_endfb71.xml"
    chain_file = "/Users/richardanderson/Downloads/chain_endfb71_pwr.xml"
    operator = openmc.deplete.CoupledOperator(model, chain_file)
    integrator = openmc.deplete.PredictorIntegrator(operator, burn_times, power)
    
    integrator.integrate()

if __name__ == "__main__":
    burn_times = [1, 10, 30, 50]  # Days
    run_burnup_simulation(enrichment=5.0, power=10e6, burn_times=burn_times, output_dir="burnup_results")

