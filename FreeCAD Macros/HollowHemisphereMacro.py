import Part
import FreeCAD as App
import FreeCADGui as Gui
import os

# Parameters
R_outer = 10      # Outer radius
thickness = 0.1     # Wall thickness
filename = "HollowHemisphere.step"
doc = App.newDocument("HollowHemisphere")

# Create full sphere and cut in half
outer_sphere = Part.makeSphere(R_outer)
half_outer = outer_sphere.cut(Part.makeBox(2*R_outer, 2*R_outer, R_outer, App.Vector(-R_outer, -R_outer, -R_outer)))

# Create inner sphere (for shell subtraction)
R_inner = R_outer - thickness
inner_sphere = Part.makeSphere(R_inner)
half_inner = inner_sphere.cut(Part.makeBox(2*R_inner, 2*R_inner, R_inner + thickness, App.Vector(-R_inner, -R_inner, -R_inner)))

# Subtract to make hollow shell
hollow_shell = half_outer.cut(half_inner)

# Add shape to document
part_obj = doc.addObject("Part::Feature", "HollowHemisphere")
part_obj.Shape = hollow_shell

# Count and print number of faces
print("Number of faces in hollow hemisphere:", len(hollow_shell.Faces))