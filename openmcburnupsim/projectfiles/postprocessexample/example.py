#!/usr/bin/env python
# coding: utf-8

# # Post Processing
# This notebook demonstrates some basic post-processing tasks that can be performed with the Python API, such as plotting a 2D mesh tally and plotting neutron source sites from an eigenvalue calculation. The problem we will use is a simple reflected pin-cell.

# In[1]:


from IPython.display import Image
import numpy as np
import matplotlib.pyplot as plt
import openmc


# ## Generate Input Files

# First we need to define materials that will be used in the problem. We'll create three materials for the fuel, water, and cladding of the fuel pin.

# In[2]:


# 1.6 enriched fuel
fuel = openmc.Material(name='1.6% Fuel')
fuel.set_density('g/cm3', 10.31341)
fuel.add_nuclide('U235', 3.7503e-4)
fuel.add_nuclide('U238', 2.2625e-2)
fuel.add_nuclide('O16', 4.6007e-2)

# borated water
water = openmc.Material(name='Borated Water')
water.set_density('g/cm3', 0.740582)
water.add_nuclide('H1', 4.9457e-2)
water.add_nuclide('O16', 2.4732e-2)
water.add_nuclide('B10', 8.0042e-6)

# zircaloy
zircaloy = openmc.Material(name='Zircaloy')
zircaloy.set_density('g/cm3', 6.55)
zircaloy.add_nuclide('Zr90', 7.2758e-3)


# With our three materials, we can now create a materials file object that can be exported to an actual XML file.

# In[3]:


# Instantiate a Materials collection
materials = openmc.Materials([fuel, water, zircaloy])

# Export to "materials.xml"
materials.export_to_xml()


# Now let's move on to the geometry. Our problem will have three regions for the fuel, the clad, and the surrounding coolant. The first step is to create the bounding surfaces -- in this case two cylinders and six reflective planes.

# In[4]:


# Create cylinders for the fuel and clad
fuel_outer_radius = openmc.ZCylinder(x0=0.0, y0=0.0, r=0.39218)
clad_outer_radius = openmc.ZCylinder(x0=0.0, y0=0.0, r=0.45720)

# Create boundary planes to surround the geometry
min_x = openmc.XPlane(x0=-0.63, boundary_type='reflective')
max_x = openmc.XPlane(x0=+0.63, boundary_type='reflective')
min_y = openmc.YPlane(y0=-0.63, boundary_type='reflective')
max_y = openmc.YPlane(y0=+0.63, boundary_type='reflective')
min_z = openmc.ZPlane(z0=-0.63, boundary_type='reflective')
max_z = openmc.ZPlane(z0=+0.63, boundary_type='reflective')


# With the surfaces defined, we can now create cells that are defined by intersections of half-spaces created by the surfaces.

# In[5]:


# Create a Universe to encapsulate a fuel pin
pin_cell_universe = openmc.Universe(name='1.6% Fuel Pin')

# Create fuel Cell
fuel_cell = openmc.Cell(name='1.6% Fuel')
fuel_cell.fill = fuel
fuel_cell.region = -fuel_outer_radius
pin_cell_universe.add_cell(fuel_cell)

# Create a clad Cell
clad_cell = openmc.Cell(name='1.6% Clad')
clad_cell.fill = zircaloy
clad_cell.region = +fuel_outer_radius & -clad_outer_radius
pin_cell_universe.add_cell(clad_cell)

# Create a moderator Cell
moderator_cell = openmc.Cell(name='1.6% Moderator')
moderator_cell.fill = water
moderator_cell.region = +clad_outer_radius
pin_cell_universe.add_cell(moderator_cell)


# OpenMC requires that there is a "root" universe. Let us create a root cell that is filled by the pin cell universe and then assign it to the root universe.

# In[6]:


# Create root Cell
root_cell = openmc.Cell(name='root cell')
root_cell.fill = pin_cell_universe

# Add boundary planes
root_cell.region = +min_x & -max_x & +min_y & -max_y & +min_z & -max_z

# Create root Universe
root_universe = openmc.Universe(universe_id=0, name='root universe')
root_universe.add_cell(root_cell)


# We now must create a geometry that is assigned a root universe, put the geometry into a geometry file, and export it to XML.

# In[7]:


# Create Geometry and set root Universe
geometry = openmc.Geometry(root_universe)


# In[8]:


# Export to "geometry.xml"
geometry.export_to_xml()


# With the geometry and materials finished, we now just need to define simulation parameters. In this case, we will use 10 inactive batches and 90 active batches each with 5000 particles.

# In[9]:


# OpenMC simulation parameters
settings = openmc.Settings()
settings.batches = 100
settings.inactive = 10
settings.particles = 5000

# Create an initial uniform spatial source distribution over fissionable zones
bounds = [-0.63, -0.63, -0.63, 0.63, 0.63, 0.63]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
settings.source = openmc.IndependentSource(space=uniform_dist)

# Export to "settings.xml"
settings.export_to_xml()


# Let us also create a plot file that we can use to verify that our pin cell geometry was created successfully.

# In[10]:


plot = openmc.Plot.from_geometry(geometry)
plot.pixels = (250, 250)
plot.to_ipython_image()


# As we can see from the plot, we have a nice pin cell with fuel, cladding, and water! Before we run our simulation, we need to tell the code what we want to tally. The following code shows how to create a 2D mesh tally.

# In[11]:


# Instantiate an empty Tallies object
tallies = openmc.Tallies()


# In[12]:


# Create mesh which will be used for tally
mesh = openmc.RegularMesh()
mesh.dimension = [100, 100]
mesh.lower_left = [-0.63, -0.63]
mesh.upper_right = [0.63, 0.63]

# Create mesh filter for tally
mesh_filter = openmc.MeshFilter(mesh)

# Create mesh tally to score flux and fission rate
tally = openmc.Tally(name='flux')
tally.filters = [mesh_filter]
tally.scores = ['flux', 'fission']
tallies.append(tally)


# In[13]:


# Export to "tallies.xml"
tallies.export_to_xml()


# Now we a have a complete set of inputs, so we can go ahead and run our simulation.

# In[14]:


# Run OpenMC!
openmc.run()


# ## Tally Data Processing

# Our simulation ran successfully and created a statepoint file with all the tally data in it. We begin our analysis here loading the statepoint file and 'reading' the results. By default, data from the statepoint file is only read into memory when it is requested. This helps keep the memory use to a minimum even when a statepoint file may be huge.

# In[15]:


# Load the statepoint file
sp = openmc.StatePoint('statepoint.100.h5')


# Next we need to get the tally, which can be done with the ``StatePoint.get_tally(...)`` method.

# In[16]:


tally = sp.get_tally(scores=['flux'])
print(tally)


# The statepoint file actually stores the sum and sum-of-squares for each tally bin from which the mean and variance can be calculated as described [here](../methods/tallies.rst#variance). The sum and sum-of-squares can be accessed using the ``sum`` and ``sum_sq`` properties:

# In[17]:


tally.sum


# However, the mean and standard deviation of the mean are usually what you are more interested in. The Tally class also has properties ``mean`` and ``std_dev`` which automatically calculate these statistics on-the-fly.

# In[18]:


print(tally.mean.shape)
(tally.mean, tally.std_dev)


# The tally data has three dimensions: one for filter combinations, one for nuclides, and one for scores. We see that there are 10000 filter combinations (corresponding to the 100 x 100 mesh bins), a single nuclide (since none was specified), and two scores. If we only want to look at a single score, we can use the ``get_slice(...)`` method as follows.

# In[19]:


flux = tally.get_slice(scores=['flux'])
fission = tally.get_slice(scores=['fission'])
print(flux)


# To get the bins into a form that we can plot, we can simply change the shape of the array since it is a numpy array.

# In[20]:


flux.std_dev.shape = (100, 100)
flux.mean.shape = (100, 100)
fission.std_dev.shape = (100, 100)
fission.mean.shape = (100, 100)


# In[21]:


fig = plt.subplot(121)
fig.imshow(flux.mean)
fig2 = plt.subplot(122)
fig2.imshow(fission.mean)


# Now let's say we want to look at the distribution of relative errors of our tally bins for flux. First we create a new variable called ``relative_error`` and set it to the ratio of the standard deviation and the mean, being careful not to divide by zero in case some bins were never scored to.

# In[22]:


# Determine relative error
relative_error = np.zeros_like(flux.std_dev)
nonzero = flux.mean > 0
relative_error[nonzero] = flux.std_dev[nonzero] / flux.mean[nonzero]

# distribution of relative errors
ret = plt.hist(relative_error[nonzero], bins=50)


# ## Source Sites

# Source sites can be accessed from the ``source`` property. As shown below, the source sites are represented as a numpy array with a structured datatype.

# In[23]:


sp.source


# If we want, say, only the energies from the source sites, we can simply index the source array with the name of the field:

# In[24]:


sp.source['E']


# Now, we can look at things like the energy distribution of source sites. Note that we don't directly use the ``matplotlib.pyplot.hist`` method since our binning is logarithmic.

# In[25]:


# Create log-spaced energy bins from 1 keV to 10 MeV
energy_bins = np.logspace(3,7)

# Calculate pdf for source energies
probability, bin_edges = np.histogram(sp.source['E'], energy_bins, density=True)

# Make sure integrating the PDF gives us unity
print(sum(probability*np.diff(energy_bins)))

# Plot source energy PDF
plt.semilogx(energy_bins[:-1], probability*np.diff(energy_bins), drawstyle='steps')
plt.xlabel('Energy (eV)')
plt.ylabel('Probability/eV')


# Let's also look at the spatial distribution of the sites. To make the plot a little more interesting, we can also include the direction of the particle emitted from the source and color each source by the logarithm of its energy.

# In[26]:


plt.quiver(sp.source['r']['x'], sp.source['r']['y'],
           sp.source['u']['x'], sp.source['u']['y'],
           np.log(sp.source['E']), cmap='jet', scale=20.0)
plt.colorbar()
plt.xlim((-0.5,0.5))
plt.ylim((-0.5,0.5))

plt.show()
# In[27]:


# Close the statepoint file as a matter of best practice
sp.close()
