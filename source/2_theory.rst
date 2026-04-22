.. _theory:

Theory of FTLE and LCS
======================

Finite-Time Lyapunov Exponents (FTLE) and Lagrangian Coherent Structures (LCS) provide a
finite-time description of transport, stirring, and material organization in fluid flows.
They are especially useful when the velocity field is known only on a bounded time interval, as is
typically the case for numerical simulations, laboratory measurements, and geophysical products.
In this setting, one is not primarily interested in instantaneous streamlines, but in the geometry
of trajectories over a prescribed interval of time. Classical invariant objects, such as stable and
unstable manifolds, still provide the conceptual background, but the practical questions are
finite-time and data-driven [HallerYuan2000]_ [Shadden2005]_ [Haller2015]_.

Mathematical Framework
----------------------

Consider an unsteady velocity field :math:`\mathbf{u}(\mathbf{x},t)` on a domain
:math:`U \subset \mathbb{R}^{n}` and the trajectory equation

.. math::

   \dot{\mathbf{x}} = \mathbf{u}(\mathbf{x},t), \qquad
   \mathbf{x}(t_0) = \mathbf{x}_0.

The associated flow map

.. math::

   \varphi_{t_0}^{t_0+T} : \mathbf{x}_0 \mapsto \mathbf{x}(t_0+T; t_0, \mathbf{x}_0)

transports an initial condition :math:`\mathbf{x}_0` from time :math:`t_0` to
:math:`t_0 + T`. The deformation of an infinitesimal perturbation
:math:`\delta \mathbf{x}_0` is described to leading order by the deformation gradient
:math:`\nabla \varphi_{t_0}^{t_0+T}` and the right Cauchy-Green strain tensor

.. math::

   \mathbf{C}_{t_0}^{t_0+T}(\mathbf{x}_0)
   =
   \left[\nabla \varphi_{t_0}^{t_0+T}(\mathbf{x}_0)\right]^{\mathsf{T}}
   \nabla \varphi_{t_0}^{t_0+T}(\mathbf{x}_0).

This tensor is symmetric and positive definite whenever the flow map is locally invertible.
Let its eigenpairs satisfy

.. math::

   \mathbf{C}_{t_0}^{t_0+T}\,\boldsymbol{\xi}_i
   =
   \lambda_i\,\boldsymbol{\xi}_i,
   \qquad
   0 < \lambda_1 \le \cdots \le \lambda_n.

Then :math:`\lambda_n = \lambda_{\max}` gives the largest finite-time stretching factor, while
:math:`\boldsymbol{\xi}_n` indicates the direction of maximal stretching at the initial time.
Accordingly, the FTLE field is defined by

.. math::

   \sigma_{t_0}^{T}(\mathbf{x}_0)
   =
   \frac{1}{|T|}
   \ln \sqrt{\lambda_{\max}\!\left(\mathbf{C}_{t_0}^{t_0+T}(\mathbf{x}_0)\right)}
   =
   \frac{1}{2|T|}
   \ln \lambda_{\max}\!\left(\mathbf{C}_{t_0}^{t_0+T}(\mathbf{x}_0)\right).

Hence FTLE measures the average exponential rate at which two initially nearby particles can
separate over the finite interval :math:`[t_0,t_0+T]`. For :math:`T>0`, the field reveals strongest
forward-time repulsion; for :math:`T<0`, it reveals strongest backward-time repulsion, which is
equivalently strongest forward-time attraction. Because it is derived from the spectrum of the
Cauchy-Green tensor, FTLE is an objective scalar diagnostic under time-dependent Euclidean
changes of frame [Shadden2005]_ [Haller2015]_.

In modern usage, an LCS is not merely a region of visually coherent trajectories, but a
codimension-one material set that organizes nearby tracer motion. In two-dimensional flows these
sets are material curves; in three-dimensional flows they are material surfaces. FTLE supplies a
useful stretching diagnostic, whereas LCS theory seeks the material geometry responsible for that
stretching [Haller2001]_ [Haller2011]_ [Haller2015]_.

.. _steady:

Steady LCS
----------

For an autonomous velocity field :math:`\dot{\mathbf{x}}=\mathbf{u}(\mathbf{x})`, the flow map forms
a one-parameter dynamical system, and the relevant organizing structures are the classical
invariant sets of nonlinear dynamics. If :math:`\mathbf{x}^{\ast}` is a hyperbolic equilibrium, its
stable and unstable manifolds are

.. math::

   W^{s}(\mathbf{x}^{\ast})
   =
   \left\{
      \mathbf{x}_0 \in U :
      \varphi_{t_0}^{t}(\mathbf{x}_0) \to \mathbf{x}^{\ast}
      \text{ as } t \to +\infty
   \right\},

.. math::

   W^{u}(\mathbf{x}^{\ast})
   =
   \left\{
      \mathbf{x}_0 \in U :
      \varphi_{t_0}^{t}(\mathbf{x}_0) \to \mathbf{x}^{\ast}
      \text{ as } t \to -\infty
   \right\}.

These manifolds are exact material barriers. In two dimensions they act as separatrices that divide
the flow into dynamically distinct regions; in three dimensions they generalize to invariant curves
and surfaces depending on the local saddle structure [HallerYuan2000]_ [Haller2001]_.

From the FTLE viewpoint, steady hyperbolic manifolds appear as sharp ridges when the integration
time :math:`|T|` is sufficiently long, because initial conditions on opposite sides of the manifold
experience substantially different fates. In this special setting, FTLE does not introduce a new
type of coherence so much as provide a finite-time visualization of an already existing invariant
geometry. Closed streamlines and invariant tori, by contrast, are associated with relatively weak
net stretching and therefore do not produce the same hyperbolic ridge signature.

This steady picture remains the correct conceptual limit for simple time-periodic or quasi-periodic
flows as well: when recurrent motion is genuinely present for all times, LCS theory asymptotically
connects with classical stable and unstable manifolds, KAM-type barriers, and other invariant
objects from dynamical systems theory [HallerYuan2000]_ [Shadden2005]_ [Haller2015]_.

.. _unsteady:

Unsteady LCS
------------

Realistic transport problems are rarely autonomous or recurrent. In aperiodic flows, and in data
sets available only on a finite interval, asymptotic notions such as stable manifolds, unstable
manifolds, or periodic orbits are generally not available as exact objects. The central question is
therefore finite-time: which material curves or surfaces organize separation, attraction, folding,
entrainment, or jet-like transport over the observation window
:math:`[t_0,t_0+T]`? [Shadden2005]_ [Haller2015]_.

The classical FTLE-based answer is to identify LCS candidates as ridges of the FTLE field.
In the influential formulation of Shadden, Lekien, and Marsden, these ridges act as finite-time
mixing templates: forward FTLE ridges approximate repelling structures, whereas backward FTLE
ridges approximate attracting structures [Shadden2005]_. This viewpoint is practically powerful
because it reduces complex trajectory behavior to a scalar field derived from the flow map.
It also explains why FTLE is widely used in oceanography, atmospheric transport, and
experimental fluid mechanics.

At the same time, a ridge of FTLE is only a diagnostic signature, not a complete material
definition of coherence. Strong ridges often mark important transport barriers, but ridge extraction
depends on the time interval, the spatial resolution, and the particular ridge criterion. Moreover,
high FTLE values can arise from strong shear without identifying a uniquely most repelling
material surface. For this reason, later work placed LCS theory on a stricter variational basis
[Haller2011]_ [Farazmand2012]_ [Haller2015]_.

In the variational theory of hyperbolic LCS, a repelling LCS is defined as a material surface whose
finite-time normal repulsion is locally maximal among nearby material surfaces; an attracting LCS
is obtained as the backward-time counterpart [Haller2011]_. This formulation links admissible
LCSs directly to the eigenvalues and eigenvectors of the Cauchy-Green tensor. In two-dimensional
flows, repelling and attracting LCSs can be constructed from special tensor lines of that field,
which explains why the eigenstructure of :math:`\mathbf{C}_{t_0}^{t_0+T}` is more fundamental than
the FTLE scalar alone [Farazmand2012]_.

The broader finite-time theory also distinguishes different transport mechanisms. Hyperbolic LCSs
govern strongest attraction and repulsion; elliptic LCSs bound vortex-like regions that resist
filamentation; parabolic LCSs act as generalized jet cores with minimal cross-stream transport
[Haller2015]_. FTLE is most naturally tied to the hyperbolic family, because it measures
exponential stretching. It is therefore an excellent first diagnostic for separation and attraction, but
it is not by itself a complete theory for all coherent transport barriers.

One further point is essential in unsteady problems: the time interval is part of the definition of
the object. Changing :math:`t_0` or :math:`T` changes the flow map, the Cauchy-Green tensor, and
hence the resulting FTLE field and LCSs. A sliding-window analysis over a long record is useful,
but each window defines a distinct finite-time dynamical system; structures extracted from
different windows should not automatically be interpreted as a single invariant object that simply
moves in time [Haller2015]_. This dependence on the chosen interval is not a defect, but a direct
reflection of the finite-time nature of observed transport.

In summary, FTLE provides an objective scalar measure of finite-time stretching, while LCS theory
seeks the material skeleton that gives this stretching geometric meaning. The FTLE-ridge picture
offers an accessible and often informative first approximation, whereas modern variational theory
clarifies when those stretching features genuinely act as repelling, attracting, vortical, or
jet-defining transport barriers [Shadden2005]_ [Haller2011]_ [Farazmand2012]_ [Haller2015]_.
