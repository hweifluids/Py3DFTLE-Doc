PyFTLE3D: 3D LCS Research Toolkit with GUI
===================================================

.. note::

   **The toolkit is undergoing an internal testing, and currently not publicly released. Who would like to use it in advance to release for non-profit use as early birds, please feel free contact me. The planned release time is 2025.06 - 2025.07 with high-freedom opensource license type, while the author still holds the right to change the release time. **
   
   **Terms marked with ``I.P.`` flag are under development. If they will be released with the toolkit at the same time is without my control due to my limited time and effort.** 
   
   **- Huanxia with Regards.**


PyFTLE3D is a Python-based research toolkit for computing and visualizing Lagrangian coherent structures (LCS) in three-dimensional fluid flows and other complex systems by computing Finite Time Lyapunov Exponent (FTLE) field.
Coming with GUI and integrated runtime, it is designed to be user-friendly and accessible, allowing peers to easily analyze and visualize LCS in their data, with no need for researchers to have any programming knowledge, just like using postprocessing software such as ParaView or Tecplot.
It provides multiple OOTB numerical methods, with various computation density, accuracy, and focus preferences.
Certainly, it is also suitable for those with programming experience who want to customize their analysis and define more functions.

The computation is accelerated by CUDA (Python Numba implement), achieving over 30x speedup compared to the CPU version.
The author currently lacks the capability to develop with OpenCL, so for now only NVIDIA GPUs are supported for acceleration.
As most computers have larger RAM than GPU memory, the toolkit also supports CPU mode for large datasets on CPU systems.
Although, the new features would firstly be added to the GPU version, and the CPU version would be updated with a month-scale delay.

For the start, it is strongly suggested to begin with reading the review by Prof. George Haller in regard to LCS `here <https://www.georgehaller.com/reprints/annurev-fluid-010313-141322.pdf>`__, published on Annual Review of Fluid Mechanics, instead of directly use this toolkit without any theoretical foundation.
Some of the theories can be found in this section, but it is still more encouraged to learn LCS theories and numerical methods by reading the technical papers.




Table of Contents
-----------------

.. toctree::
   :maxdepth: 2

   index
   1_requirements
   2_theory
   3_numerical
   4_input
   5_gui
   6_inputdeck
   7_command
   8_faqs
   9_references


.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
