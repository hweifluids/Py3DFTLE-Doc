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

The first-order explicit Euler scheme ``Euler`` advances the position by sampling the velocity at the beginning of the time step:

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

**Trilinear**

The ``trilinear`` interpolation method is the fastest among all the methods given in ``PyFTLE3D``, which is a low-order implement.
The continuous velocity field is reconstructed by trilinear interpolation of the component maps ``u``, ``v``, ``w`` that live on the eight vertices of a Cartesian cell:

.. math::

   i   = \left\lfloor\frac{x-x_0}{\Delta x}\right\rfloor,\;
   j   = \left\lfloor\frac{y-y_0}{\Delta y}\right\rfloor,\;
   k   = \left\lfloor\frac{z-z_0}{\Delta z}\right\rfloor,\\
   \tau_x = \frac{x-x_0}{\Delta x}-i,\;
   \tau_y = \frac{y-y_0}{\Delta y}-j,\;
   \tau_z = \frac{z-z_0}{\Delta z}-k.

.. math::

   \mathbf{u}=
   \sum_{d_x,d_y,d_z\in\{0,1\}}
   (1-d_x+(-1)^{d_x}\tau_x)
   (1-d_y+(-1)^{d_y}\tau_y)
   (1-d_z+(-1)^{d_z}\tau_z)
   \mathbf{u}_{\,i+d_x,\,j+d_y,\,k+d_z}.




**Tricubic Catmull-Rom**

The Tricubic Catmull-Rom interpolation ``tricubic`` used here is a separable three-dimensional cubic spline based on the one-dimensional Catmul-Rom spline (parameter ``a=-0.5``) for velocity fields, which ensures :math:`C^{1}` continuity.
The process is given as follows.

The one-dimensional Catmull-Rom interpolation reconstructs a :math:`C^{1}`-continuous approximation of velocity at an arbitrary location :math:`x = i + t`, where :math:`t \in [0,1)` and :math:`i = \lfloor x\rfloor` on a uniform grid with :math:`\Delta x = 1`. 
A four-point stencil is used:

.. math::

   \{\,u_{i-1},\,u_{i},\,u_{i+1},\,u_{i+2}\}.

Define coefficients for :math:`t \in [0,1)`:

.. math::

   a_{0} &= -\tfrac{1}{2}\,u_{i-1} + \tfrac{3}{2}\,u_{i} - \tfrac{3}{2}\,u_{i+1} + \tfrac{1}{2}\,u_{i+2},  \\[6pt]
   a_{1} &= u_{i-1} - 2.5\,u_{i} + 2.0\,u_{i+1} - 0.5\,u_{i+2},  \\[6pt]
   a_{2} &= -0.5\,u_{i-1} + 0.5\,u_{i+1},  \\[6pt]
   a_{3} &= u_{i}.

The one-dimensional interpolant is:

.. math::

   u_{\mathrm{CR}}(i + t) 
   = ((a_{0}\,t + a_{1})\,t + a_{2})\,t + a_{3}.

For three-dimensional interpolation, let the target velocity location be

.. math::

   (f_{x},\,f_{y},\,f_{z}), 
   \quad
   i = \lfloor f_{x}\rfloor,\;
   j = \lfloor f_{y}\rfloor,\; 
   k = \lfloor f_{z}\rfloor,
   \quad
   t_{x} = f_{x} - i,\;
   t_{y} = f_{y} - j,\;
   t_{z} = f_{z} - k.

At each :math:`z = k + \Delta` (where :math:`\Delta \in \{-1,0,1,2\}`), perform bicubic interpolation (first in ``x``, then in ``y``). For each :math:`y = j + \ell` (where :math:`\ell \in \{-1,0,1,2\}`), compute:

.. math::

   M_{\,j+\ell}(k+\Delta)
   = \mathrm{CR}_{1}\bigl(
     u_{\,i-1,\,j+\ell,\,k+\Delta},\,
     u_{\,i,\,j+\ell,\,k+\Delta},\,
     u_{\,i+1,\,j+\ell,\,k+\Delta},\,
     u_{\,i+2,\,j+\ell,\,k+\Delta}
     ;\,t_{x}\bigr).

Then combine along ``y``:

.. math::

   B(k+\Delta)
   = \mathrm{CR}_{1}\bigl(
     M_{\,j-1}(k+\Delta),\,
     M_{\,j}(k+\Delta),\,
     M_{\,j+1}(k+\Delta),\,
     M_{\,j+2}(k+\Delta)
     ;\,t_{y}\bigr).

Finally, interpolate along the third direction ``z``:

.. math::

   u(f_{x},\,f_{y},\,f_{z})
   = \mathrm{CR}_{1}\bigl(
     B(k-1),\,B(k),\,B(k+1),\,B(k+2)
     ;\,t_{z}\bigr).

Here :math:`\mathrm{CR}_{1}(⋯; t)` denotes the one-dimensional Catmull-Rom spline defined above.



**Tricubic by F. Lekien** ``I.P.``

Marked as ``tricubicFL``, this variation of tricubic interpolator is still under development by the author.

**Hermite**

The ``hermite`` 

**WENO**

The weighted essentially non-oscillatory ``WENO`` used here is a fifth-order WENO reconstruction (WENO-5). It is suggested to be used in research with intermittent capture need, e.g., high-speed flows and shock capture.
It shows relatively poor performance in general cases, and comsuming more wall time. 
The method originates from [Jiang1996]_ and expanded to three-dimensional computation, and [Shu2009]_ gave a review on the WEMO method.
The process is given as follows.

The WENO-5 method reconstructs a non-oscillatory, fifth-order-accurate approximation of a function value at an arbitrary location :math:`x = x_{i+1/2} + t\,\Delta x`, where :math:`t \in [0,1)` and :math:`x_{i+1/2} = x_i + \tfrac{1}{2}\,\Delta x` on a uniform grid with :math:`\Delta x = 1`. 
A five-point stencil ``f_{i-2}``, ``f_{i-1}``, ``f_i``, ``f_{i+1}``, ``f_{i+2}`` is used.

.. math::

   \{\,f_{i-2}, f_{i-1}, f_{i}, f_{i+1}, f_{i+2}\}.

Define three overlapping three-point stencils:

.. math::

   S_{0} = \{f_{i-2}, f_{i-1}, f_{i}\}, \quad
   S_{1} = \{f_{i-1}, f_{i}, f_{i+1}\}, \quad
   S_{2} = \{f_{i}, f_{i+1}, f_{i+2}\}.

On each stencil :math:`S_{\ell}` (:math:`\ell = 0,1,2`), construct a quadratic polynomial 

.. math::

   p_{\ell}(t) = C_{\ell,0} + C_{\ell,1}\,t + C_{\ell,2}\,t^{2}, \quad \ell = 0,1,2,

that interpolates the three values in that stencil at :math:`x = x_{i+1/2} + t\,\Delta x`.

The coefficients are chosen so that each :math:`p_{\ell}(t)` matches :math:`f` at the three stencil points. 

For 
:math:`S_{0} = \{f_{i-2}, f_{i-1}, f_{i}\}`:

.. math::

   C_{0,0} = \frac{2\,f_{i-2} - 7\,f_{i-1} + 11\,f_{i}}{6},\\
   C_{0,1} = \frac{-f_{i-2} + 5\,f_{i-1} - 4\,f_{i} + f_{i+1}}{2},\\
   C_{0,2} = \frac{f_{i-2} - 2\,f_{i-1} + f_{i}}{2}.

For :math:`S_{1} = \{f_{i-1}, f_{i}, f_{i+1}\}`:

.. math::

   C_{1,0} = \frac{-f_{i-1} + 5\,f_{i} + 2\,f_{i+1}}{6},\\
   C_{1,1} = \frac{f_{i-1} - f_{i+1}}{2},\\
   C_{1,2} = \frac{f_{i-1} - 2\,f_{i} + f_{i+1}}{2}.

For :math:`S_{2} = \{f_{i}, f_{i+1}, f_{i+2}\}`:

.. math::

   C_{2,0} = \frac{2\,f_{i} + 5\,f_{i+1} - f_{i+2}}{6},\\
   C_{2,1} = \frac{-f_{i} + 4\,f_{i+1} - 3\,f_{i+2}}{2},\\
   C_{2,2} = \frac{f_{i} - 2\,f_{i+1} + f_{i+2}}{2}.

Once :math:`p_{0}(t)`, :math:`p_{1}(t)`, and :math:`p_{2}(t)` are defined, compute the Jiang–Shu smoothness indicators :math:`\beta_{\ell}` for each stencil:

.. math::

   \beta_{0} = 13\,\bigl(f_{i-2} - 2\,f_{i-1} + f_{i}\bigr)^{2} + 3\,\bigl(f_{i-2} - 4\,f_{i-1} + 3\,f_{i}\bigr)^{2},\\
   \beta_{1} = 13\,\bigl(f_{i-1} - 2\,f_{i} + f_{i+1}\bigr)^{2} + 3\,\bigl(f_{i-1} - f_{i+1}\bigr)^{2},\\
   \beta_{2} = 13\,\bigl(f_{i} - 2\,f_{i+1} + f_{i+2}\bigr)^{2} + 3\,\bigl(3\,f_{i} - 4\,f_{i+1} + f_{i+2}\bigr)^{2}.

Fixed linear weights are :math:`(d_{0}, d_{1}, d_{2}) = (0.1,\,0.6,\,0.3)`. Introduce :math:`\varepsilon = 10^{-6}` and define unnormalized weights:

.. math::

   \tilde{\alpha}_{\ell} = \frac{d_{\ell}}{(\varepsilon + \beta_{\ell})^{2}}, \quad \ell = 0,1,2.

Normalize to obtain nonlinear weights :math:`\omega_{\ell}`:

.. math::

   \omega_{\ell} = \frac{\tilde{\alpha}_{\ell}}{\tilde{\alpha}_{0} + \tilde{\alpha}_{1} + \tilde{\alpha}_{2}}, \quad \sum_{\ell=0}^{2} \omega_{\ell} = 1.

Finally, reconstruct at :math:`x = x_{i+1/2} + t\,\Delta x` by combining:

.. math::

   f_{\mathrm{WENO5}}(x_{i+1/2} + t\,\Delta x) = \omega_{0}\,p_{0}(t) + \omega_{1}\,p_{1}(t) + \omega_{2}\,p_{2}(t).




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


.. _numdyn:

Windowing for Dynamic LCS
----------------------------------------------



.. _numcompare:

Computational Density and Comparison
----------------------------------------------

.. _numtips:

General Tips
~~~~~~~~~~~~~~~~


As for your reference, and configured as defaults, the *Berkeley LCS Tutorials* used ``RK4`` for advection.
The velocity fields were interpolated with ``tricubic-FL`` by them, originating from [Lekien2005]_, which has higher performance by solving a 64×64 linear system using the function values, gradients, and mixed partial derivatives at its eight corners, which is in future development plan for ``Py3DFTLE`` with high priority.
Although not detailed, ``grad_order=2`` was employed by them from the equation, supposing the mesh is sufficiently refined.

Please always notice that, although providing much better numerical precision and looks cool in papers, high-order methods could be resource-consuming, even several hundred times.