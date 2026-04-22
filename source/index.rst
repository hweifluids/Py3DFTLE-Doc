Streamcenter+: Native CUDA/Qt Toolkit for FTLE and LCS
========================================================

.. note::

   The toolkit is still undergoing internal testing and is currently not publicly released.
   Those who would like to use it early for non-profit work are welcome to contact the author.
   The previously announced release window was 2025.06-2025.07, but the schedule may still change.

   Terms marked with ``I.P.`` are under development and may not be released at the same time as
   the main toolkit.

   Huanxia with regards.


Streamcenter+ is a GPU-accelerated toolbox for computing and visualizing Lagrangian coherent
structures (LCS) in two- and three-dimensional unsteady flows through Finite-Time Lyapunov
Exponent (FTLE) fields.
The current mainline is centered on a native CUDA solver and a Qt/VTK desktop GUI, combining
structured VTK time-sequence input, ``.par``-driven execution, and integrated result inspection
in one workflow.

The previous Python implementation is preserved as a legacy archive for reference and migration,
but it is no longer the supported mainline path.
For the mainline solver, NVIDIA GPUs are currently the supported acceleration target.
Multiple advection, interpolation, wall-treatment, and dynamic-window options are provided so the
user can balance computational density, smoothness, and numerical cost for a given research task.

For the start, it is strongly suggested to begin with reading the review by Prof. George Haller in regard to LCS `here <https://www.georgehaller.com/reprints/annurev-fluid-010313-141322.pdf>`__, published on Annual Review of Fluid Mechanics, instead of directly use this toolkit without any theoretical foundation.
Some of the theories can be found in this section, but it is still more encouraged to learn LCS theories and numerical methods by reading the technical papers.




Table of Contents
-----------------

.. toctree::
   :maxdepth: 2

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
