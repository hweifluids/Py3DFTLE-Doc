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

Time Integration
~~~~~~~~~~~~~~~~~~

Consider the initial‐value problem for passive tracer advection in a continuous velocity field

.. math::

   \frac{d\mathbf{x}}{dt} = \sigma\,\mathbf{u}(\mathbf{x},t)\,,  
   \mathbf{x}(t_n)=\mathbf{x}_n\,,  

where :math:`\sigma = \pm1` selects forward or backward integration.

**Explicit Euler Method**

The first‐order explicit Euler scheme advances the position by sampling the velocity at the beginning of the time step:

.. math::

   \mathbf{u}_n = \mathbf{u}(\mathbf{x}_n,t_n),\\
   \mathbf{x}_{n+1} = \mathbf{x}_n + \sigma\,\Delta t\,\mathbf{u}_n.

This method incurs a global error of order :math:`O(\Delta t)` and requires only one velocity evaluation per step.

**Second‐Order Runge–Kutta (Heun’s Method)**

Heun’s method attains second‐order accuracy by combining predictor and corrector slopes:

.. math::

   k_1 = \sigma\,\mathbf{u}(\mathbf{x}_n,t_n),\\
   \mathbf{x}^* = \mathbf{x}_n + \Delta t\,k_1,\\
   k_2 = \sigma\,\mathbf{u}(\mathbf{x}^*,t_n + \Delta t),\\
   \mathbf{x}_{n+1} = \mathbf{x}_n + \tfrac{\Delta t}{2}\,(k_1 + k_2).

This scheme yields a global error of order :math:`O(\Delta t^2)` with two velocity evaluations per step.

**Classical Fourth‐Order Runge–Kutta (RK4)**

The classical RK4 method achieves fourth‐order accuracy via four slope evaluations at intermediate points:

.. math::

   k_1 = \mathbf{u}(\mathbf{x}_n,t_n),\\
   k_2 = \mathbf{u}\!\bigl(\mathbf{x}_n + \tfrac{\Delta t}{2}k_1,\;t_n + \tfrac{\Delta t}{2}\bigr),\\
   k_3 = \mathbf{u}\!\bigl(\mathbf{x}_n + \tfrac{\Delta t}{2}k_2,\;t_n + \tfrac{\Delta t}{2}\bigr),\\
   k_4 = \mathbf{u}(\mathbf{x}_n + \Delta t\,k_3,\;t_n + \Delta t),\\
   \mathbf{x}_{n+1} = \mathbf{x}_n + \tfrac{\Delta t}{6}\,(k_1 + 2k_2 + 2k_3 + k_4).

This yields a global error of order :math:`O(\Delta t^4)` with four velocity evaluations per step.

**Sixth‐Order Runge–Kutta (RK6)**

The six‐stage scheme uses non‐uniform weights to attain sixth‐order accuracy:

.. math::

   k_1 = \mathbf{u}(\mathbf{x}_n,t_n),\\
   k_2 = \mathbf{u}\!\bigl(\mathbf{x}_n + \tfrac{\Delta t}{3}k_1,\;t_n + \tfrac{\Delta t}{3}\bigr),\\
   k_3 = \mathbf{u}\!\bigl(\mathbf{x}_n + \Delta t(\tfrac{1}{6}k_1 + \tfrac{1}{6}k_2),\;t_n + \tfrac{\Delta t}{3}\bigr),\\
   k_4 = \mathbf{u}\!\bigl(\mathbf{x}_n + \Delta t(\tfrac{1}{8}k_1 + \tfrac{3}{8}k_3),\;t_n + \tfrac{\Delta t}{2}\bigr),\\
   k_5 = \mathbf{u}\!\bigl(\mathbf{x}_n + \Delta t(\tfrac{1}{2}k_1 - \tfrac{3}{2}k_3 + 2k_4),\;t_n + \tfrac{2\Delta t}{3}\bigr),\\
   k_6 = \mathbf{u}\!\bigl(\mathbf{x}_n + \Delta t(-\tfrac{3}{2}k_1 + 2k_2 - \tfrac{1}{2}k_3 + k_4),\;t_n + \Delta t\bigr),\\
   \mathbf{x}_{n+1} = \mathbf{x}_n + \Delta t\bigl(\tfrac{1}{20}k_1 + \tfrac{1}{4}k_4 + \tfrac{1}{5}k_5 + \tfrac{1}{2}k_6\bigr).

This scheme incurs a global error of order :math:`O(\Delta t^6)` with six velocity evaluations.  

All methods assume a continuous velocity interpolation (e.g., tricubic) to supply :math:`\mathbf{u}` at arbitrary particle positions and times.

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