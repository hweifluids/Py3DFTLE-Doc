.. _numerical:

Numerical Methods
===================

.. _advc:

Lagrangian Advection
-----------------------



.. _intp:
Velocity Interpolation
~~~~~~~~~~~~~~~~~~~~~~~

.. _wall:

Wall Treatment
~~~~~~~~~~~~~~~

.. _marching:

Particle Marching
~~~~~~~~~~~~~~~~~~



.. _ftlefinal:

FTLE Computation
-------------------








.. _grad:
Gradient Discretization
~~~~~~~~~~~~~~~~~~~~~~~~~




.. _eigen:
Eigenvalue Solver
~~~~~~~~~~~~~~~~~~~~




.. _numcompare:

Computational Density and Comparison
==========================================

.. _numtips:

General Tips
~~~~~~~~~~~~~~~~


As for your reference, and configured as defaults, the *Berkeley LCS Tutorials* used ``RK4`` for advection.
The velocity fields were interpolated with ``tricubic-FL`` by them, originating from [Lekien2005]_, which has higher performance by solving a 64×64 linear system using the function values, gradients, and mixed partial derivatives at its eight corners, which is in future development plan for ``Py3DFTLE`` with high priority.
Although not detailed, ``grad_order=2`` was employed by them from the equation, supposing the mesh is sufficiently refined.

Please always notice that, although providing much better numerical precision and looks cool in papers, high-order methods could be resource-consuming, even several hundred times.