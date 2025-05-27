PyFTLE3D: Research Toolkit for 3D Lagrangian coherent structures with GUI
===================================================

**Notice: The toolkit is undergoing an internal testing, and currently not publicly released. Who would like to use it in advance to release for non-profit use as early birds, please feel free contact me. The planned release time is 2025.06 - 2025.07 with high-freedom opensource license type, while the author still holds the right to change the release time. - Huanxia with Regards**

PyFTLE3D is a Python-based research toolkit for computing and visualizing
Lagrangian coherent structures (LCS) in three-dimensional fluid flows and other complex systems by computing Finite Time Lyapunov Exponent.
Coming with GUI and integrated environment, it is designed to be user-friendly and accessible, allowing peers to easily analyze and visualize LCS in their data
with multiple OOTB numerical methods, with various computation density, accuracy, and focus preferences,
with no need for researchers to have any programming knowledge.
Certainly, it is also suitable for those with programming experience who want to customize their analysis.

The computation is accelerated by CUDA (Python Numba), the NVIDIA treasure, achieving over 30x speedup compared to the CPU version.
The author currently lacks the capability to develop with OpenCL, so for now only NVIDIA GPUs are supported for acceleration.
As most of computers have larger RAM than GPU memory, the toolkit also supports CPU mode for large datasets on CPU systems.
Although, the new features would firstly be added to the GPU version, and the CPU version would be updated with month-scale delay.



This guide is intended to help new users get started with using nekRS, as well as serve as a
reference for more advanced users. Because the :term:`Nek5000` code is somewhat of a predecessor to
nekRS, some aspects of the current nekRS design are selected to enable faster translation of
:term:`Nek5000` input files into nekRS input files. Throughout this documentation, all such
:term:`Nek5000`-oriented settings will be referred to as "legacy" settings. Because these
:term:`Nek5000`-oriented settings require proficiency in Fortran, structured text formats,
and several additional input files, all new users are encouraged to adopt the nekRS-based problem setup.

We recommend working through this user guide in the order below. At the very least, please
read :ref:`The nekRS Input Files <input>` page before reading the :ref:`FAQs <detailed>`
page, as some necessary concepts are introduced in this order.

.. note::

   This documentation is a work in progress, and will undergo big changes as more
   features are added to nekRS. Please open issues to track any missing information
   at the github repository `here <https://github.com/aprilnovak/nekRSDoc>`__.

Table of Contents
-----------------

.. toctree::
   :maxdepth: 2

   requirements
   theory
   input_files
   commonly_used_variables
   compiling
   tutorials
   detailed_usage
   plugins
   just_in_time_compilation
   doxygen
   glossary
   references

.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
