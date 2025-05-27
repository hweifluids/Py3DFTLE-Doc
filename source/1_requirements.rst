.. _requirements:

Requirements and Quickstart
========================

Basic Runtime
-----------------

This page summarizes the system requirements and dependencies of NekRS.

1. Linux or Mac OS X (Microsoft WSL and Windows is not supported)
2. GNU/oneAPI/NVHPC/ROCm compilers (C++17/C99 compatible)
3. MPI version 3.1 or later
4. [CMake](https://cmake.org/) version 3.21 or later


Add-on Libs
-----------------
1. ParaView
ParaView is a powerful open-source visualization tool that can be used to visualize the results of 3D FTLE computations. We also integrated the ParaView entrance into our GUI, so you can directly open ParaView from it.

Quickstart
-----------------
An computation example coming with a set of large eddy simulation (LES) data is provided.
The dataset can be downloaded from `here <https://noting.noting>`, and the case descriptions can be referred to paper [1]_.
