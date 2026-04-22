.. _input:

Input Data Files
=====================

This page describes what constitutes a qualified velocity-field input for the current native
solver. In both GUI and command-line workflows, the solver reads a folder of time-ordered
legacy ``.vtk`` files given by ``vtk_folder`` in :ref:`inputdeck`.
Each file represents one time frame of the same flow field, and the whole folder is interpreted as
one temporal sequence.


Supported File Type
-------------------

The current native solver scans only files whose extension is ``.vtk``.
These files are loaded through VTK's legacy-data reader. As a consequence:

- legacy VTK ``.vtk`` is the supported input container for the solver;
- ASCII and binary legacy files are both acceptable as long as VTK can parse them;
- XML VTK files such as ``.vts`` or ``.vti`` are not used as solver input in the current mainline,
  even though they may be viewable elsewhere.


Folder-Based Time Sequence
--------------------------

The solver treats the input as a folder-level sequence rather than a single file:

- only regular ``.vtk`` files in ``vtk_folder`` are considered;
- subfolders are not scanned recursively;
- files are sorted by the last integer appearing in the filename;
- this sorted order defines the frame index used by ``T_range``.

For this reason, filenames such as
``frame_00000.vtk``, ``frame_00001.vtk``, ``frame_00002.vtk`` are strongly recommended.
Missing frame numbers are allowed, but the ordering still follows the extracted integer tags.
If a filename contains more than one integer group, the last one is the one used for sorting.


Structured Spatial Layout
-------------------------

The native solver assumes that every frame lies on one structured tensor-product lattice in
:math:`x`, :math:`y`, and :math:`z`.
In practice, the most reliable qualified inputs are legacy VTK files whose dataset is either
``RECTILINEAR_GRID`` or ``STRUCTURED_GRID``.

For a sequence to be qualified:

- all frames must describe the same spatial grid;
- the point count must remain identical across all frames;
- the coordinate axes must be consistent from frame to frame;
- each point must be mappable to one unique lattice index :math:`(i,j,k)`.

For ``RECTILINEAR_GRID``, the solver takes the axes directly from
``X_COORDINATES``, ``Y_COORDINATES``, and ``Z_COORDINATES``.
For other structured datasets, the solver reconstructs the axes from the point coordinates.
Therefore, even when ``STRUCTURED_GRID`` is used, the point coordinates should still lie on
separable :math:`x`, :math:`y`, :math:`z` axes.
General unstructured meshes are not qualified input for the present solver.


Velocity Arrays
---------------

Velocity must be stored in ``POINT_DATA``.
Arrays provided only in ``CELL_DATA`` are not sufficient for the current solver.

The parser currently accepts any of the following point-data conventions:

- three scalar arrays named ``u``, ``v``, ``w``;
- three scalar arrays named ``U``, ``V``, ``W``;
- three scalar arrays named ``Ux``, ``Uy``, ``Uz``;
- three scalar arrays named ``Vx``, ``Vy``, ``Vz``;
- two scalar arrays named ``u`` and ``v``;
- two scalar arrays named ``U`` and ``V``;
- two scalar arrays named ``Ux`` and ``Uy``;
- two scalar arrays named ``Vx`` and ``Vy``;
- one vector array named ``U``, ``V``, ``Velocity``, ``velocity``, ``Vel``, ``vel``, ``data``,
  or ``Data`` with at least two components.

If only two components are present, the solver fills the third component with zero and therefore
treats the input as quasi-two-dimensional.
All frames in the sequence should use the same array naming convention.


Minimal Qualified Skeleton
--------------------------

A minimal qualified legacy-VTK frame may look conceptually like this:

.. code-block:: text

   # vtk DataFile Version 5.1
   velocity frame
   BINARY
   DATASET RECTILINEAR_GRID
   DIMENSIONS Nx Ny Nz
   X_COORDINATES Nx double
   ...
   Y_COORDINATES Ny double
   ...
   Z_COORDINATES Nz double
   ...
   POINT_DATA Nx*Ny*Nz
   VECTORS Velocity double
   ...

An equally acceptable scalar-array variant is:

.. code-block:: text

   POINT_DATA Nx*Ny*Nz
   SCALARS u double
   LOOKUP_TABLE default
   ...
   SCALARS v double
   LOOKUP_TABLE default
   ...
   SCALARS w double
   LOOKUP_TABLE default
   ...

The exact numeric precision is not prescribed here, but the data must be readable by VTK and the
point-data array lengths must match the total number of grid points.


Two-Dimensional and Thin-Slab Data
----------------------------------

The current native solver can also work with quasi-two-dimensional inputs.
In that case:

- the :math:`x` and :math:`y` directions must still contain at least two grid points;
- the :math:`z` direction may be a single layer or a very thin slab;
- the third velocity component may be omitted if the input follows one of the accepted
  two-component conventions above.

This is useful for planar velocity fields embedded in a three-dimensional file format.


Practical Checklist
-------------------

Before launching the solver, it is recommended to confirm the following:

- ``vtk_folder`` points to a directory, not to a single file;
- the directory contains only the intended ``.vtk`` time frames, or at least no unrelated
  ``.vtk`` files whose filename numbers would disturb the ordering;
- all frames share one structured grid and one consistent point-data layout;
- velocity is stored in ``POINT_DATA``;
- the sorted sequence is long enough for the requested ``T_range``.

If ``auto_wall = 1``, the solver uses the outer coordinates of the input grid as the seeding bounds.
If ``auto_wall = 0``, the wall ranges in the ``.par`` file must instead define the computational
domain explicitly.
