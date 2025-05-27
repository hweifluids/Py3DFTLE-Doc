.. _requirements:

Requirements and Quickstart
========================

Basic Environment
-----------------

This page summarizes the system requirements and dependencies of NekRS.

1. `numpy <https://numpy.org>`__ version 1.21 or later  
2. `matplotlib <https://matplotlib.org>`__ version 3.4 or later  
3. `pyvista <https://pyvista.org>`__ version 0.32 or later  
4. `scipy <https://scipy.org>`__ version 1.7 or later  
5. `pyevtk <https://github.com/paulo-herrera/PyEVTK>`__ version 1.4 or later  
6. `numba <https://numba.pydata.org>`__ version 0.54 or later  
7. `tqdm <https://tqdm.github.io>`__ version 4.62 or later  
8. `cupy-cuda12x <https://cupy.dev>`__  
9. `pyvistaqt <https://github.com/pyvista/pyvistaqt>`__ version 0.11.2 or later  
10. `PyQt5 <https://riverbankcomputing.com/software/pyqt/intro>`__ version 5.15.11 or later  



Add-on Libs
-----------------
1. **ParaView**

ParaView is a powerful beloved open-source visualization platform based on vtk that can be used to visualize the results of 3D FTLE computations. We also integrated the ParaView entrance into our GUI, so you can directly open ParaView from it. The ParaView installation-free package (v6.0) can be downloaded from `official <https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v6.0&type=binary&os=Windows&downloadFile=ParaView-6.0.0-RC1-MPI-Windows-Python3.12-msvc2017-AMD64.zip>`__. 

Quickstart
-----------------
An computation example coming with a set of our in-house downscaled Direct Navier-Stokes Simulation (DNS) data is provided, including :math:
The dataset can be downloaded from `here <https://noting.noting>`, and the case descriptions can be referred to paper [1]_.
