.. _requirements:

Requirements and Quickstart
========================

Hardware and Hardware Requirements
----------------------------------

The minimiums are listed as follows:

* Minimium 4GB RAM and 10GB free disk.
* Not sensitive to platform. Windows 11 and Ubuntu 22.04 are fully tested.

GPU acceleration is suggested for large dataset computation. Therefore, the requirements for GPU acceleration is listed as follows:

* NVIDIA GPUs, with a minimum `compute capability <https://developer.nvidia.com/cuda-gpus>`__ of 5.0, i.e., later than GTX 750.
* It is suggested that, make ``Nx`` × ``Ny`` × ``Nz`` × ``Nt`` less than 200,000,000 Lagrangian particles on each Gigabyte of GPU memory for double-side computation under default numerical methods for one still. Generally half it under high-order computations.
For example, you can execute on a 600×300×300 mesh for c.a. 80 time steps for a still.
* The GPU requirements for dynamic LCS are undergoing experiments by the author.


Basic Environment
-----------------

This page summarizes the system requirements and dependencies of ``PyFTLE3D``.

Firstly, the following packages are supposed to be installed on your computer manually in advance.

1. `Python <https://www.python.org/>`__ version 3.13. Later version than 3.8 till 3.13 should work properly, but was not fully tested. The official download could be found at `here <https://www.python.org/downloads/release/python-3130/?featured_on=pythonbytes>`__. Please select Python installer according to your system framework.
2. `CUDA Toolkit <https://developer.nvidia.com/cuda-toolkit>`__ version 12.9 (`Download <https://developer.nvidia.com/cuda-toolkit-archive>`__). All versions crossing 12.x.x should theoretically work, but only 12.9 was fully tested. Ealier versions could work as well, but it requires a different version of `cupy`, and could cause unexpected performance decay and errors. See `cupy` documentation for more information.
3. `pip <https://pypi.org/project/pip/>`__ newest version, which is used for installing further dependencies. It can be installed by:

.. code-block::

  python -m ensurepip --upgrade

Generally, the following dependencies can be installed via ``pip``, the Python package manager.

1. `numpy <https://numpy.org>`__ version 1.21 or later  
2. `matplotlib <https://matplotlib.org>`__ version 3.4 or later  
3. `pyvista <https://pyvista.org>`__ version 0.32 or later  
4. `scipy <https://scipy.org>`__ version 1.7 or later  
5. `pyevtk <https://github.com/paulo-herrera/PyEVTK>`__ version 1.4 or later  
6. `numba <https://numba.pydata.org>`__ version 0.54 or later  
7. `tqdm <https://tqdm.github.io>`__ version 4.62 or later  
8. `cupy-cuda12x <https://pypi.org/project/cupy-cuda12x/>`__ version 13.4.1 or later
9. `pyvistaqt <https://github.com/pyvista/pyvistaqt>`__ version 0.11.2 or later  
10. `PyQt5 <https://riverbankcomputing.com/software/pyqt/intro>`__ version 5.15.11 or later  

You can install them after installing ``pip`` by:

.. code-block::

  pip install -r ./requirements.txt

When it does not work, consider if your current dir is incorrect. It should be run from the project root of ``Py3DFTLE`` when undergoing configurations and under :ref:`command-line <command>` mode. Please feel free about such thing under GUI mode.


Add-on Libs
-----------------

**1. ParaView** is a powerful beloved open-source visualization platform based on *vtk* that can be used to visualize the results of 3D FTLE computations. We also integrated the ParaView entrance into our GUI, so you can directly open ParaView from it.
The ParaView installation-free package (v6.0) can be downloaded from `official <https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v6.0&type=binary&os=Windows&downloadFile=ParaView-6.0.0-RC1-MPI-Windows-Python3.12-msvc2017-AMD64.zip>`__. 
For already-installed ParaView, set environment variable ``PARAVIEW_PATH`` pointing to the root dir of ParaView to enable integration from GUI of ``Py3DFTLE``.

Quickstart
-----------------

An computation example coming with a set of our in-house downscaled Direct Navier-Stokes Simulation (DNS) data is provided, including :math:
The dataset can be downloaded from `here <https://noting.noting>`__, and the case descriptions can be referred to paper ``I.P.``.
