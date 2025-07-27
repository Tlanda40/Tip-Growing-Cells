# Tip-Growing-Cells
This repository is a work in progress! It simulates the cell wall expansion
of a tip growing cell using a finite-element framework in MATLAB. The most
recent iteration uses a single-domain geometry and the femodel toolbox. Previous
iterations have used multi-domain geometries and other frameworks, such as the
structural pdemodel and general pdemodel toolboxes. Files used in the most
recent iterations are directly exposed in Tip-Growing-Cells, whereas older
versions, references, and other tools are stored in subfolders in case they are
needed. When files in subfolders are no longer needed, they are removed. See the
Table of Contents for more information.

# TABLE OF CONTENTS
# Rico's Code
This model is a hardcoded 2D description of cell wall expansion at a tip,
created by Rico. It is a reference for this project.

# Reference
These files, created by Dylan, are simple geometries undergoing expansion. They
were used as references in the creation of more complex pdemodels (and also to
teach me how FE works).

# Structural PDEMODEL
These files use the static-structural pdemodel toolbox to create basic,
high-level models. They were used to test certain geometric procedures and more
complex Geometries created natively and/or imported from Free CAD or Blender.
These files also contain an initial implementation of turgor pressure
renormalization. Extensibility is uniform.

# General PDEMODEL
These files use the 3D pdemodel toolbox to create low-level, high-control
models. Material property tensors are explicitly defined. Unfortunately, the
solver breaks frequently.

# Multidomain FEMODEL
These files use the femodel toolbox to create a relatively high-level
simulation. The geometry is multidomain, enabling greater control and variable
extensibility, but with significant runtimes (approximately one hour in total).
Certain matrices have been cached in the subfolder to make the runtime less than
a second in testing.

# FreeCAD Macros
These files are FreeCad macros (and their CAD realizations) used to create the
initial geometry.

# STEP Geometries
These files are initial geometries created in FreeCad.
