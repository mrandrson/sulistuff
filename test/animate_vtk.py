#!/usr/bin/env pvpython
import glob
from paraview.simple import *

files = sorted(glob.glob("post_vtk/dump*.vtk"))

reader = OpenDataFile(files)
reader.UseAllTimes = True
reader.UpdatePipeline()

view = GetActiveViewOrCreate('RenderView')
display = Show(reader, view)

glyph = Glyph(reader, GlyphType='Sphere')
glyph.ScaleMode = 'Scale by Array'
glyph.ScaleArray = ['POINT_DATA', 'radius']
glyph.ScaleFactor = 0.01  # adjust to your scale
Show(glyph, view)

Render()
for t in reader.TimestepValues:
    view.ViewTime = t
    Render()
    # SaveScreenshot(f"frame_{int(t):04d}.png", view)

