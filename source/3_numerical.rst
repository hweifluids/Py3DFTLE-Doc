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

Consider the initial-value problem for particle advection in a velocity field

.. math::

   \frac{d\mathbf{x}}{dt} = \sigma\,\mathbf{u}(\mathbf{x},t),
   \quad \mathbf{x}(t_n) = \mathbf{x}_n,

where :math:`\sigma = \pm 1` indicates forward/backward advection.

**Explicit Euler Method**

Given :math:`\mathbf{x}_n` at time :math:`t_n`, the first-order (Euler) update is

.. math::

   \mathbf{x}_{n+1} = \mathbf{x}_n + \Delta t\,f(\mathbf{x}_n,t_n)
                     = \mathbf{x}_n + \sigma\,\Delta t\,\mathbf{u}(\mathbf{x}_n,t_n)

Implementation steps:

  1. Evaluate velocity:  
     .. math::  
        \mathbf{u}_n = \mathbf{u}(\mathbf{x}_n, t_n)  
  2. Advance position:  
     .. math::  
        \mathbf{x}_{n+1} = \mathbf{x}_n + \sigma\,\Delta t\,\mathbf{u}_n  

This method is simple but incurs :math:`O(\Delta t)` local truncation error.

**Second-Order Runge–Kutta (Heun’s Method)**

A two-stage explicit scheme with :math:`O(\Delta t^2)` accuracy:

.. math::

   k_1 &= f(\mathbf{x}_n, t_n),\\
   \mathbf{x}^* &= \mathbf{x}_n + \Delta t\,k_1,\\
   k_2 &= f(\mathbf{x}^*, t_n + \Delta t),\\
   \mathbf{x}_{n+1} &= \mathbf{x}_n + \tfrac{\Delta t}{2}\,(k_1 + k_2).

Implementation steps:

  1. Compute :math:`k_1 = \sigma\,\mathbf{u}(\mathbf{x}_n, t_n)`.  
  2. Predict step:  
     .. math::  
        \mathbf{x}^* = \mathbf{x}_n + \sigma\,\Delta t\,k_1  
  3. Compute :math:`k_2 = \sigma\,\mathbf{u}(\mathbf{x}^*, t_n + \Delta t)`.  
  4. Update:  
     .. math::  
        \mathbf{x}_{n+1} = \mathbf{x}_n + \tfrac{\sigma\,\Delta t}{2}\,(k_1 + k_2)

**Classical Fourth-Order Runge–Kutta (RK4)**

A four-stage scheme with :math:`O(\Delta t^4)` accuracy:

.. math::

   k_1 &= f(\mathbf{x}_n, t_n),\\
   k_2 &= f\!\bigl(\mathbf{x}_n + \tfrac{\Delta t}{2}k_1,\;t_n + \tfrac{\Delta t}{2}\bigr),\\
   k_3 &= f\!\bigl(\mathbf{x}_n + \tfrac{\Delta t}{2}k_2,\;t_n + \tfrac{\Delta t}{2}\bigr),\\
   k_4 &= f(\mathbf{x}_n + \Delta t\,k_3,\;t_n + \Delta t),\\
   \mathbf{x}_{n+1} &= \mathbf{x}_n + \tfrac{\Delta t}{6}\,(k_1 + 2k_2 + 2k_3 + k_4).

Implementation steps:

  1. Compute :math:`k_1` at :math:`\mathbf{x}_n`.  
  2. Compute :math:`k_2` at :math:`\mathbf{x}_n + \tfrac{\Delta t}{2}k_1`.  
  3. Compute :math:`k_3` at :math:`\mathbf{x}_n + \tfrac{\Delta t}{2}k_2`.  
  4. Compute :math:`k_4` at :math:`\mathbf{x}_n + \Delta t\,k_3`.  
  5. Combine:  
     .. math::  
        \mathbf{x}_{n+1} = \mathbf{x}_n + \tfrac{\Delta t}{6}\,(k_1 + 2k_2 + 2k_3 + k_4)

**Sixth-Order Runge–Kutta (RK6)**

A six-stage explicit method with :math:`O(\Delta t^6)` accuracy. Define:

.. math::

   k_1 &= f(\mathbf{x}_n, t_n),\\
   k_2 &= f\!\bigl(\mathbf{x}_n + \tfrac{\Delta t}{3}k_1,\;t_n + \tfrac{\Delta t}{3}\bigr),\\
   k_3 &= f\!\Bigl(\mathbf{x}_n + \Delta t\bigl(\tfrac{1}{6}k_1 + \tfrac{1}{6}k_2\bigr),\;t_n + \tfrac{\Delta t}{3}\Bigr),\\
   k_4 &= f\!\Bigl(\mathbf{x}_n + \Delta t\bigl(\tfrac{1}{8}k_1 + \tfrac{3}{8}k_3\bigr),\;t_n + \tfrac{\Delta t}{2}\Bigr),\\
   k_5 &= f\!\Bigl(\mathbf{x}_n + \Delta t\bigl(\tfrac{1}{2}k_1 - \tfrac{3}{2}k_3 + 2k_4\bigr),\;t_n + \tfrac{2\Delta t}{3}\Bigr),\\
   k_6 &= f\!\Bigl(\mathbf{x}_n + \Delta t\bigl(-\tfrac{3}{2}k_1 + 2k_2 - \tfrac{1}{2}k_3 + k_4\bigr),\;t_n + \Delta t\Bigr),\\
   \mathbf{x}_{n+1} &= \mathbf{x}_n + \Delta t\Bigl(\tfrac{1}{20}k_1 + \tfrac{1}{4}k_4 + \tfrac{1}{5}k_5 + \tfrac{1}{2}k_6\Bigr).

Implementation steps:

  1. Compute each :math:`k_i = \sigma\,\mathbf{u}(\cdot)` at its intermediate point.  
  2. Form the weighted sum:  
     .. math::  
        \mathbf{x}_{n+1} = \mathbf{x}_n + \Delta t\Bigl(\tfrac{1}{20}k_1 + \tfrac{1}{4}k_4 + \tfrac{1}{5}k_5 + \tfrac{1}{2}k_6\Bigr)

All methods assume a continuous, differentiable velocity field via tricubic interpolation; replacing :math:`f` by the chosen sampler affects only boundary‐condition treatment.


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