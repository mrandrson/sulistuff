import openmc
import shutil

fuel = openmc.Material(name='fuel')
fuel.add_nuclide('U235', 1.0)
fuel.set_density('g/cm3', 10.0)

water = openmc.Material(name='water')
water.add_nuclide('H1', 2.0)
water.add_nuclide('O16', 1.0)
water.set_density('g/cm3', 1.0)

mats = openmc.Materials([fuel, water])
mats.export_to_xml()

r_pin = openmc.ZCylinder(r=0.25)
fuel_cell = openmc.Cell(fill=fuel, region=-r_pin)
water_cell = openmc.Cell(fill=water, region=+r_pin)
pin_universe = openmc.Universe(cells=[fuel_cell, water_cell])

outer_radius = 4.0
outer_surface = openmc.ZCylinder(r=outer_radius, boundary_type='vacuum')
main_cell = openmc.Cell(fill=pin_universe, region=-outer_surface)

geom = openmc.Geometry([main_cell])
geom.export_to_xml()

bounds = geom.bounding_box
width_x = abs(bounds[1][0] - bounds[0][0])
width_y = abs(bounds[1][1] - bounds[0][1])

pixels_per_unit = 150  
pixels_x = int(width_x * pixels_per_unit)
pixels_y = int(width_y * pixels_per_unit)

plot = openmc.Plot()
plot.width = (width_x, width_y)
plot.origin = (0.0, 0.0, 0.0)
plot.color_by = "material"
plot.pixels = (pixels_x, pixels_y)

plots = openmc.Plots([plot])
plots.export_to_xml()

openmc.plot_geometry()

shutil.move("plot_1.png", "/projectfiles/example.png")

print(f"Saved plot with width: {width_x}, height: {width_y}, resolution: {pixels_x}x{pixels_y}")

