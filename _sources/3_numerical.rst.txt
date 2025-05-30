.. _numerical:

Numerical Methods
===================

.. _advc:

Lagrangian Advection
-----------------------


In the numerical computation of the ``FTLE``, we first compute the *flow map* :math:`\varphi_{t_n}^{t_{n+1}}(x_n)`, which maps the initial seed point :math:`x_n` at time :math:`t_n` to time :math:`t_{n+1}`.
To obtain this map, one must numerically integrate the underlying dynamical system, which is described by the ODE:

.. math::

   \frac{d\mathbf{x}}{dt} = \sigma\,\mathbf{u}(\mathbf{x},t)\,,  
   \mathbf{x}(t_n)=\mathbf{x}_n\,,  

where :math:`\sigma = \pm1` selects forward or backward integration.
During the time integration process, the algorithm frequently queries the flow velocity vector :math:`\mathbf{u}(\mathbf{x},t)` at specific locations and moments with very high precision requirements.
However, since the data grid is inherently spatially discretized, high-order interpolation methods are required to keep numerical stability and obtain a physically meaningful flow map.
Furthermore, when querying the velocity, special wall treatment must be applied at the boundaries to avoid value discontinuities and to represent certain real physical conditions.



.. _marching:

Time Integration
~~~~~~~~~~~~~~~~~~


**Explicit Euler Method**

The first-order explicit Euler scheme advances the position by sampling the velocity at the beginning of the time step:

.. math::

   \mathbf{u}_n = \mathbf{u}(\mathbf{x}_n,t_n),\\
   \mathbf{x}_{n+1} = \mathbf{x}_n + \sigma\,\Delta t\,\mathbf{u}_n.

This method incurs a global error of order :math:`O(\Delta t)` and requires only one velocity evaluation per step, therefore has high computational speed.

**Runge-Kutta Method**

Proposed by Carl Runge and Martin Kutta around 1900, Runge-Kutta methods constitute a widely used family of algorithms for the numerical integration of ODEs.

In an explicit :math:`s`-stage Runge-Kutta scheme for this initial-value problem, the solution is advanced over a time step :math:`\Delta t` as follows.
First, compute the intermediate stage vectors:

.. math::

   \mathbf{k}_i = 
   \mathbf{u} \Bigl(   \mathbf{x}_n + \sigma\,\Delta t \sum_{j=1}^{i-1} a_{ij}\,\mathbf{k}_j,\,
   t_n + c_i \Delta t \Bigr),
   \qquad i = 1, 2, \dots, s,

and then update the solution:

.. math::

   \mathbf{x}_{n+1} =
   \mathbf{x}_n + \sigma\,\Delta t \sum_{i=1}^s b_i\,\mathbf{k}_i.

Here, the boldface stage variables :math:`\mathbf{k}_i` represent intermediate slope estimates.


**Second-Order Runge-Kutta (RK2, Heun's)**

Heun's ``RK2`` method attains second-order accuracy by combining predictor and corrector slopes:

.. math::

   k_1 = \sigma\,\mathbf{u}(\mathbf{x}_n,t_n),\\
   \mathbf{x}^* = \mathbf{x}_n + \Delta t\,k_1,\\
   k_2 = \sigma\,\mathbf{u}(\mathbf{x}^*,t_n + \Delta t),\\
   \mathbf{x}_{n+1} = \mathbf{x}_n + \tfrac{\Delta t}{2}\,(k_1 + k_2).

This scheme yields a global error of order :math:`O(\Delta t^2)` with two velocity evaluations per step.

**Classical Fourth-Order Runge-Kutta (RK4)**

The classical ``RK4`` method achieves fourth-order accuracy via four slope evaluations at intermediate points:

.. math::

   k_1 = \mathbf{u}(\mathbf{x}_n,t_n),\\
   k_2 = \mathbf{u}\!\bigl(\mathbf{x}_n + \tfrac{\Delta t}{2}k_1,\;t_n + \tfrac{\Delta t}{2}\bigr),\\
   k_3 = \mathbf{u}\!\bigl(\mathbf{x}_n + \tfrac{\Delta t}{2}k_2,\;t_n + \tfrac{\Delta t}{2}\bigr),\\
   k_4 = \mathbf{u}(\mathbf{x}_n + \Delta t\,k_3,\;t_n + \Delta t),\\
   \mathbf{x}_{n+1} = \mathbf{x}_n + \tfrac{\Delta t}{6}\,(k_1 + 2k_2 + 2k_3 + k_4).

This yields a global error of order :math:`O(\Delta t^4)` with four velocity evaluations per step.

**Sixth-Order Runge-Kutta (RK6)**

The seven-stage scheme ``ERK6(7)`` uses non-uniform weights to attain global :math:`O(\Delta t^6)` accuracy, originating from [Butcher1964]_.
As for the coefficients for ``RK6`` are more complex to write into equations, the Butcher table is given as follows.

+-----------------------------+-------------------------------+----------------------------+---------------------------------+------------------------------+-----------------------------+--------------------------+-------------------------+
| :math:`c_i`                 | :math:`a_{i1}`                | :math:`a_{i2}`             | :math:`a_{i3}`                  | :math:`a_{i4}`               | :math:`a_{i5}`              | :math:`a_{i6}`           | :math:`a_{i7}`          |
+=============================+===============================+============================+=================================+==============================+=============================+==========================+=========================+
| :math:`0`                   |                               |                            |                                 |                              |                             |                          |                         |
+-----------------------------+-------------------------------+----------------------------+---------------------------------+------------------------------+-----------------------------+--------------------------+-------------------------+
| :math:`(5\mp\sqrt{5})/10`   | :math:`(5\mp\sqrt{5})/10`     |                            |                                 |                              |                             |                          |                         |
+-----------------------------+-------------------------------+----------------------------+---------------------------------+------------------------------+-----------------------------+--------------------------+-------------------------+
| :math:`(5\pm\sqrt{5})/10`   | :math:`\mp\sqrt{5}/10`        | :math:`(5\pm2\sqrt{5})/10` |                                 |                              |                             |                          |                         |
+-----------------------------+-------------------------------+----------------------------+---------------------------------+------------------------------+-----------------------------+--------------------------+-------------------------+
| :math:`(5\mp\sqrt{5})/10`   | :math:`(-15\pm7\sqrt{5})/20`  | :math:`(-1\pm\sqrt{5})/4`  | :math:`(15\mp7\sqrt{5})/10`     |                              |                             |                          |                         |
+-----------------------------+-------------------------------+----------------------------+---------------------------------+------------------------------+-----------------------------+--------------------------+-------------------------+
| :math:`(5\pm\sqrt{5})/10`   | :math:`(5\mp\sqrt{5})/60`     | :math:`0`                  | :math:`1/6`                     | :math:`(15\pm7\sqrt{5})/60`  |                             |                          |                         |
+-----------------------------+-------------------------------+----------------------------+---------------------------------+------------------------------+-----------------------------+--------------------------+-------------------------+
| :math:`(5\mp\sqrt{5})/10`   | :math:`(5\pm\sqrt{5})/60`     | :math:`0`                  | :math:`(9\mp5\sqrt{5})/12`      | :math:`1/6`                  | :math:`(-5\pm3\sqrt{5})/10` |                          |                         |
+-----------------------------+-------------------------------+----------------------------+---------------------------------+------------------------------+-----------------------------+--------------------------+-------------------------+
| :math:`1`                   | :math:`1/6`                   | :math:`0`                  | :math:`(-55\pm25\sqrt{5})/12`   | :math:`(-25\mp7\sqrt{5})/12` | :math:`5\mp2\sqrt{5}`       | :math:`(5\pm\sqrt{5})/2` |                         |
+-----------------------------+-------------------------------+----------------------------+---------------------------------+------------------------------+-----------------------------+--------------------------+-------------------------+
| :math:`b_i`                 | :math:`1/12`                  | :math:`0`                  | :math:`0`                       | :math:`0`                    | :math:`5/12`                | :math:`5/12`             | :math:`1/12`            |
+-----------------------------+-------------------------------+----------------------------+---------------------------------+------------------------------+-----------------------------+--------------------------+-------------------------+

In our computation, the up symbol side is applied, in other words, ``±`` represents ``+``, taking :math:`\lambda=+\sqrt{5}`. With 15 digis are kept, the explicit Butcher table for ``RK6`` used by the author is shown in the following table.

+-------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+
|  :math:`c_i`                  |  :math:`a_{i1}`               |  :math:`a_{i2}`               |  :math:`a_{i3}`               |  :math:`a_{i4}`               |  :math:`a_{i5}`               |  :math:`a_{i6}`               |  :math:`a_{i7}`               |
+===============================+===============================+===============================+===============================+===============================+===============================+===============================+===============================+
| 0                             | 0                             | 0                             | 0                             | 0                             | 0                             | 0                             | 0                             |
+-------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+
| 0.276393202250021             | 0.276393202250021             | 0                             | 0                             | 0                             | 0                             | 0                             | 0                             |
+-------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+
| 0.723606797749979             | -0.223606797749979            | 0.947213595499958             | 0                             | 0                             | 0                             | 0                             | 0                             |
+-------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+
| 0.276393202250021             | 0.0326237921249264            | 0.309016994374947             | -0.0652475842498529           | 0                             | 0                             | 0                             | 0                             |
+-------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+
| 0.723606797749979             | 0.0460655337083368            | 0                             | 0.166666666666667             | 0.510874597374975             | 0                             | 0                             | 0                             |
+-------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+
| 0.276393202250021             | 0.12060113295833              | 0                             | -0.181694990624912            | 0.166666666666667             | 0.170820393249937             | 0                             | 0                             |
+-------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+
| 1                             | 0.166666666666667             | 0                             | 0.0751416197912285            | -3.38770632020821             | 0.52786404500042              | 3.61803398874989              | 0                             |
+-------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+
| :math:`\mathbf{b_i}`          | 0.0833333333333333            | 0                             | 0                             | 0                             | 0.416666666666667             | 0.416666666666667             | 0.0833333333333333            |
+-------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+



.. _intp:

Velocity Interpolation
~~~~~~~~~~~~~~~~~~~~~~~

.. _wall:

Wall Treatment
~~~~~~~~~~~~~~~












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