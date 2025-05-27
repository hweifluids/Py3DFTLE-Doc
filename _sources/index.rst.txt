PyFTLE3D: 3D LCS Research Toolkit with GUI
===================================================

.. note::

   **Notice: The toolkit is undergoing an internal testing, and currently not publicly released. Who would like to use it in advance to release for non-profit use as early birds, please feel free contact me. The planned release time is 2025.06 - 2025.07 with high-freedom opensource license type, while the author still holds the right to change the release time.** 
   **- Huanxia with Regards.**


PyFTLE3D is a Python-based research toolkit for computing and visualizing Lagrangian coherent structures (LCS) in three-dimensional fluid flows and other complex systems by computing Finite Time Lyapunov Exponent (FTLE) field.
Coming with GUI and integrated runtime, it is designed to be user-friendly and accessible, allowing peers to easily analyze and visualize LCS in their data, with no need for researchers to have any programming knowledge, just like using postprocessing software such as ParaView or Tecplot.
It provides multiple OOTB numerical methods, with various computation density, accuracy, and focus preferences.
Certainly, it is also suitable for those with programming experience who want to customize their analysis and define more functions.

The computation is accelerated by CUDA (Python Numba implement), achieving over 30x speedup compared to the CPU version.
The author currently lacks the capability to develop with OpenCL, so for now only NVIDIA GPUs are supported for acceleration.
As most computers have larger RAM than GPU memory, the toolkit also supports CPU mode for large datasets on CPU systems.
Although, the new features would firstly be added to the GPU version, and the CPU version would be updated with a month-scale delay.

For the start, it is strongly suggested to begin with reading the review by Prof. George Haller in regard to LCS `here<https://www.georgehaller.com/reprints/annurev-fluid-010313-141322.pdf>__`, published on Annual Review of Fluid Mechanics, instead of directly use this toolkit without any theoretical foundation.
Some of the theories can be found in this section, it is more encouraged to learn LCS theories and numerical methods by reading the technical papers.


We recommend working through this user guide in the order below. At the very least, please
read :ref:`The nekRS Input Files <input>` page before reading the :ref:`FAQs <detailed>`
page, as some necessary concepts are introduced in this order.



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
