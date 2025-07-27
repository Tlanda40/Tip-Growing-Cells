import FreeCAD as App
import Part

doc = App.ActiveDocument
doc.recompute()

# === Parameters ===
solid = doc.getObject("HollowHemisphere")
nx = 20   # Number of slices in X
ny = 20   # Number of slices in Y
nz = 10   # Number of slices in Z

# === Get Bounding Box ===
bbox = solid.Shape.BoundBox
xmin, xmax = bbox.XMin, bbox.XMax
ymin, ymax = bbox.YMin, bbox.YMax
zmin, zmax = bbox.ZMin, bbox.ZMax

dx = (xmax - xmin) / nx
dy = (ymax - ymin) / ny
dz = (zmax - zmin) / nz

cell_count = 0

# === Iterate Over 3D Grid ===
for i in range(nx):
    for j in range(ny):
        for k in range(nz):
            x0 = xmin + i * dx
            y0 = ymin + j * dy
            z0 = zmin + k * dz

            # Create slicing box for this grid cell
            box = Part.makeBox(dx, dy, dz, App.Vector(x0, y0, z0))
            intersection = solid.Shape.common(box)

            if not intersection.isNull() and intersection.Volume > 0:
                obj = doc.addObject("Part::Feature", f"Cell_{i}_{j}_{k}")
                obj.Shape = intersection
                obj.ViewObject.Visibility = False
                cell_count += 1
            else:
                print(f"Skipping empty or invalid cell: Cell_{i}_{j}_{k}")
                continue

doc.recompute()
print(f"Created {cell_count} sliced cells.")