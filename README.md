# Tip-Growing-Cells
This repository is a bunch of code that simulates the cell wall expansion
of a tip growing cell using a finite-element framework in MATLAB. The most
recent iteration uses a single-domain geometry and the femodel toolbox,
alongside tools in the C++ computational geometry library CGAL. Previous
iterations have used multi-domain geometries and other frameworks, such as the
structural pdemodel and general pdemodel toolboxes. Files used in the most
recent iterations are directly exposed in Tip-Growing-Cells, whereas older
versions, references, and other tools are stored in subfolders in case they are
needed. When files in subfolders are no longer needed, they are removed. See the
**Table of Contents** for descriptions of subfolder contents, see **Exposed Files**
for descriptions  and see **Instructions** for help with setup, installation
of dependencies, and troubleshooting. This repo is a work in progress! See
**Updates** to see my work pipeline and the latest updates.
---
### `main.mpp`
## Table of Contents
### FreeCAD Macros
These files are FreeCad macros (and their CAD realizations) used to create the
initial geometry.

### General PDEMODEL
These files use the 3D pdemodel toolbox to create low-level, high-control
models. Material property tensors are explicitly defined. Unfortunately, the
solver breaks frequently.

### Multidomain FEMODEL
These files use the femodel toolbox to create a relatively high-level
simulation. The geometry is multidomain, enabling greater control and variable
extensibility, but with significant runtimes (approximately one hour in total).
Certain matrices have been cached in the subfolder to make the runtime less than
a second in testing.

### Reference
These files, created by Dylan, are simple geometries undergoing expansion. They
were used to teach me how FE works in matlab and in general.

### Rico's Code
This model is a hardcoded 2D description of cell wall expansion at a tip,
created by Rico. It is a reference for this project. To view the model, run
`SimulateTipGrowth12280.m`.

### STEP Geometries
These (currently unused) files are initial geometries created in FreeCad.

### Structural PDEMODEL
These files use the static-structural pdemodel toolbox to create basic,
high-level models. They were used to test certain geometric procedures and more
complex Geometries created natively and/or imported from Free CAD or Blender.
These files also contain an initial implementation of turgor pressure
renormalization. Extensibility is uniform.

### cgal_geodesic
This is the C++ source directory. It contains the source code which implements
CGAL's MMP exact geodesic in `geodesic_example.cpp`, as well as the configuration
file and the build directory, which contains mesh data, compiled files, and the
executable. This is *not* an unused subfolder. It is the only folder which
contains files that are called by the files currently exposed in this repo.

### Exposed Files

