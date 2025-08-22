# Tip-Growing-Cells
This repository is a bunch of code that simulates the cell wall expansion
of a tip growing cell using a finite-element framework in MATLAB. The most
recent iteration uses a single-domain geometry and the femodel toolbox,
alongside tools in the C++ computational geometry library CGAL. Previous
iterations have used multi-domain geometries and other frameworks, such as the
structural pdemodel and general pdemodel toolboxes. Files used in the most
recent iterations are directly exposed in Tip-Growing-Cells, whereas older
versions, references, and other tools are stored in the subfolder **References**
in case they are needed. When files in subfolders are no longer needed, they are
removed. See the **Table of Contents** for descriptions of the subfolders and 
exposed files and see **Instructions** for help with setup and installation of 
dependencies. This repo is a work in progress!

# Table of Contents
### FreeCAD Macros
These files are FreeCad macros used to create the
initial geometry.

### References
#### General PDEMODEL
These files use the 3D pdemodel toolbox to create low-level, high-control
models. Material property tensors are explicitly defined. Unfortunately, the
solver breaks frequently.

#### Multidomain FEMODEL
These files use the femodel toolbox to create a relatively high-level
simulation. The geometry is multidomain, enabling greater control and variable
extensibility, but with significant runtimes (approximately one hour in total).
Certain matrices have been cached in the subfolder to make the runtime less than
a second in testing.

#### Basic References
These files, created by Dylan, are simple geometries undergoing expansion. They
were used to teach me how FE works in matlab and in general.

#### STEP Geometries
These (currently unused) files are initial geometries created in FreeCad.

#### Structural PDEMODEL
These files use the static-structural pdemodel toolbox to create basic,
high-level models. They were used to test certain geometric procedures and more
complex Geometries created natively and/or imported from Free CAD or Blender.
These files also contain an initial implementation of turgor pressure
renormalization. Extensibility is uniform.

### Rico's Code
This model is a hardcoded 2D description of cell wall expansion at a tip,
created by Rico. It is a reference for this project. To view the model, run
`SimulateTipGrowth12280.m`.

### cgal_geodesic
This is the C++ source directory. It contains the source code which implements
CGAL's MMP exact geodesic in `geodesic_example.cpp`, as well as the configuration
file and the build directory, which contains mesh data, compiled files, and the
executable. This is *not* an unused subfolder. It is the only folder which
contains files that are called by the files currently exposed in this repo. In
particular, it is called by `FEMODELVERSIONMMP.m`.

### `ExtensibilityProfileTester.m`
This file is used to plot and tune extensibility profiles and Young's Modulus
profiles.

### `FEMODELVERSIONBASIC.m`
This file is the simplest implementation of a single-domain FEMODEL with
nonconstant material properties. Is varies material properties with xy-distance
instead of geodesic distance, for simplicity, and is used as a reference for
the other exposed FEMODEL implementations. It also has a time-loop implemented, for
tip growth.

### `FEMODELVERSIONDIJKSTRA.m`
This is an implementation of FEMODEL which uses Dijkstra's algorithm to
calculate geodesic distances within the mesh. Material properties are assigned
correctly, but their values are less accurate (overestimated) than the values
in `FEMODELVERSIONMMP.m`. Because Dijkstra's algorithm is implemented
natively in MATLAB, this file is easy to run, and it is also low runtime.
Because of these useful properties, it is used for development and
troubleshooting of `FEMODELVERSIONMMP.m`.

### `FEMODELVERSIONMMP.m`
This is an implementation of FEMODEL which uses the Mitchell-Mount-Papadimitriou
exact geodesic algorithm (MMP) to calculate geodesic distances within the mesh.
It is the most advanced and accurate model in the repo. It also runs the C++
executable in `/cgal_geodesic/build`.

### `HollowHemisphere.step`
This is a .step file of the initial geometry, created in freeCAD. This geometry
is a hemispherical shell with a bottom. Its thickness is one one-hundredth of
its radius. This file is called in all exposed FEMODEL implementations, and
most of the archived files in subfolders.

### `README.md`
The user guide, written in Markdown!

### `write_off.m`
This file writes a .off file given matrices V and F (vertices and faces) of a mesh.
It is called by `FEMODELVERSIONMMP.m` in order to export the mesh into a format
acceptable to CGAL functions.

# Instructions
To ensure that files in subfolders are accessible to exposed files, it is necessary
to add them to the path. Change the path to `Tip-Growing-Cells`, then run
`addpath(genpath(pwd)))` in the MATLAB console. Then install dependencies, build,
and run. Bear in mind that the only MATLAB file that requires the C++ requirements
and setup is `FEMODELVERSIONMMP.m`, so if you just want to run other things, those
are not necessary.

## 🧱 Dependencies
To build and run this project, you'll need:

### C++ Requirements
- [CMake](https://cmake.org/download/) (for the C++ build process)
- A C++ compiler (GCC or Clang for macOS; Visual Studio for Windows)
- [CGAL](https://www.cgal.org/download/) (The Computational Geometry Algorithms Library)
- [Boost](https://www.boost.org/) — required by some CGAL functions
- [Eigen](https://eigen.tuxfamily.org/) — required by CGAL

### MATLAB Requirements
- This code has only been run on MATLAB R2025a or later. Using significantly earlier versions
may cause errors.
- Two special MATLAB toolboxes are required: "Statistics and Machine Learning Toolbox" and "Partial Differential Equation Toolbox". Other toolboxes may be required by some older archived versions.

## C++ Setup
### 🪟 Windows (with Visual Studio)

1. Install [CMake](https://cmake.org/download/) and [Visual Studio](https://visualstudio.microsoft.com/) with C++ tools.
2. Open the "x64 Native Tools Command Prompt for VS".
3. Run the following commands:

```bash
cd cgal_geodesic
mkdir build
cd build
cmake .. -G "Visual Studio 18 2025" 
cmake --build . --config Release
```

### 🍎 macOS / 🐧 Linux

1. Install CMake and a C++ compiler (Clang for macOS, or GCC)
2. Open Terminal
3. Build the code:

```bash
cd cgal_geodesic
mkdir build
cd build
cmake ..
make
```

✅ These processes create the executable inside `CGAL/build/`. Then you're good to run `FEMODELVERSIONMMP.m`!
Note that MATLAB behaves mostly the same across OSes, except when calling executables. However, the MATLAB
script auto-detects OS and chooses between `geodesic_example` (mac/Linux) or `geodesic_example.exe` (Windows),
so no changes need to be made to the MATLAB file.
