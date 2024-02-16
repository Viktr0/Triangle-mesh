# Triangle Mesh

This is a minimal framework for 3D mesh manipulation, which can be extended and used as a base for various projects, for example prototypes for fairing algorithms, or even displaying/modifying parametric surfaces, etc.

A framework for practicing 
* handling of 3D points, vectors and triangle meshes
* displaying points, sections and polygons
* coloring, texturing
* some basic algorithms (curvature calculation, smoothing, evaluating Bezier surfaces)

<p align="center">
  <img width=70% src="https://github.com/Viktr0/Triangle-mesh/assets/47856193/e6bce3e3-4423-4c0e-ac63-9ba2490ed9cd"/>
</p>

## Usage

The following hotkeys are available:
* `H`: Displays a help window that contains all the necessary information.
* `R`: Reload model
* `O`: Toggle orthographic projection
* `P`: Set plain map (no coloring)
* `M`: Set mean curvature map
* `L`: Set slicing map
     - +: Increase slicing density
     - -: Decrease slicing density
     - *: Set slicing direction to view
* `I`: Set isophote line map
* `E`: Set environment texture
* `C`: Toggle control polygon visualization
* `S`: Toggle solid (filled polygon) visualization
* `W`: Toggle wireframe visualization
* `X`: Colors the triangles whose area is "typical" of the triangle mesh, i.e. larger than the smallest 25% and smaller than the largest 25% (ie between the 1st and 3rd quartiles), and also writes the value of the median.
* `F`: Fair mesh

<p align="center">
  <img src="https://github.com/Viktr0/Triangle-mesh/assets/47856193/a21e81a9-55e9-48e1-8804-b6ab20bd908d" alt="animated" />
</p>

## Installations

Dependencies:

- Qt5 (http://qt-project.org/)
- libQGLViewer (http://www.libqglviewer.com/)
- OpenMesh (http://www.openmesh.org/)

### Linux

    qmake && make

Assumes that:

- qmake (Qt 5.x) is in the path
- QGLViewer include files are in `<standard include directory>/QGLViewer/`
- OpenMesh  include files are in `<standard include directory>/OpenMesh/`
- QGLViewer library file  is  in `<standard library directory>/`
- OpenMesh  library files are in `/usr/lib/OpenMesh/`

If any of the above is not satisfied, edit sample-framework.pro accordingly.

### Windows / Visual Studio

1. Install the Qt SDK, which should integrate itself into Visual Studio.

1. Download the source package for libQGLViewer, and put it somewhere,
   e.g. at c:\Program Files\libQGLViewer. Open Visual Studio,
   in the Qt menu select "Open Qt project file (*.pro)",
   and open the project file in the QGLViewer subdirectory
   (in my case c:\Program Files\libQGLViewer\libQGLViewer.pro).
   Compile a release version of the library.

1. Now open sample-framework.pro, and modify the following line:

        LIBQGLVIEWER_INSTALL_PATH = 'C:\Program Files\libQGLViewer'
   (using the correct path on your system, of course).

1. Download the source package for OpenMesh, and put it somewhere,
   e.g. at c:\Program Files\OpenMesh. Open the solution file
   in Visual Studio and build a release version of the core library.
   Then open sample-framework.pro, and modify the following line:

        OPENMESH_INSTALL_PATH = 'C:\Program Files\OpenMesh'
   (using the correct path on your system, of course).

1. Open Visual Studio, in the Qt menu select "Open Qt project file (*.pro)",
   and open sample-framework.pro.


1. You should be able to build the project, but it won't start. Solution:
   copy QGLViewer2.dll (found in QGLViewer\QGLViewer\release\)
   into c:\Windows\System32 or into the project's directory.

1. If you also want to build a debug version of sample-framework,
   you are still not ready! You have to build debug versions of libQGLViewer
   and OpenMesh first, then change the library names in the project properties
   dialog window (and don't forget to copy QGLViewerd2.dll to a location
   in the path).
