commercial :term:`CFD<CFD>` codes
is ``Nelements * Np``, 
the :math:`x`-coordinates  a :math:`\checkmark` in 
in the
:ref:`Creating a Mesh for Conjugate Heat Transfer <cht_mesh>` section

Plugins
=======

RANS :math:`k`-:math:`\tau` Plugin
----------------------------------

Add the Closure Properties Calculation
""""""""""""""""""""""""""""""""""""""
  1. Set the ``udf.sEqnSource`` function pointer to a function
     local to the ``.udf`` file that actually computes the source terms
  2. Add that source term function to the ``.udf``

.. code-block::
  user$ export NEKRS_HOME=$HOME/.local/nekrs
  
.. code-block:: cpp

  void user_q(nrs_t *nrs, dfloat time, occa::memory o_S, occa::memory o_FS)
  {
    RANSktau::updateSourceTerms();
  }


.. warning::

  nekRS's :math:`k`-:math:`\tau` implementation currently requires that
  the laminar dynamic viscosity and the density are constant. Therefore, you
  should not have any other material properties being set in this function
  like there were in :ref:`Setting Custom Properties <custom_properties>`.

.. _setting_ICs:
Then, be sure to add this directory to your path:

    ALE
      Arbitrary Lagrange Eulerian
ed in [Persson]_
from `github <https://github.com/Nek5000/nekRS>`__.

================== ============================ ================== =================================================
Variable Name      Size                         Device?            Meaning
================== ============================ ================== =================================================
``comm``           1                                               MPI communicator
``device``         1                                               backend device
``dim``            1                                               spatial dimension of mesh
``elementInfo``    ``Nelements``                                   phase of element (0 = fluid, 1 = solid)
``EToB``           ``Nelements * Nfaces``       :math:`\checkmark` boundary ID for each face
``N``              1                                               polynomial order for each dimension
``NboundaryFaces`` 1                                               *total* number of faces on a boundary (rank sum)
``Nelements``      1                                               number of elements
``Nfaces``         1                                               number of faces per element
``Nfp``            1                                               number of quadrature points per face
``Np``             1                                               number of quadrature points per element
``rank``           1                                               parallel process rank
``size``           1                                               size of MPI communicator
``vmapM``          ``Nelements * Nfaces * Nfp`` :math:`\checkmark` quadrature point index for faces on boundaries
``x``              ``Nelements * Np``           :math:`\checkmark` :math:`x`-coordinates of quadrature points
``y``              ``Nelements * Np``           :math:`\checkmark` :math:`y`-coordinates of quadrature points
``z``              ``Nelements * Np``           :math:`\checkmark` :math:`z`-coordinates of quadrature points
================== ============================ ================== =================================================
n :numref:`fig-walls`.
Additionally *Nek5000* can handle conjugate heat transfer problems.

.. _fig-walls:

.. figure:: figs/walls.png
    :align: center
    :figclass: align-center
    :alt: domains

    Computational domain showing respective fluid and solid subdomains, :math:`\Omega_f` and
    :math:`\Omega_s`.  The shared boundaries are denoted :math:`\partial\Omega_f=\partial\Omega_s`
    and the solid boundary which is not shared by fluid is :math:`\overline{\partial\Omega_s}`,
    while the fluid boundary not shared by solid :math:`\overline{\partial\Omega_f}`.
.. _flow_vars:

Flow Solution Fields and Simulation Settings
Flow Solution Fields and Simulation Settings
--------------------------------------------
corresponding array ``nrs->mue`` member.

================== ================================= ================== ======================================================================================================
Variable Name      Size                              Device?            Meaning
================== ================================= ================== ======================================================================================================
``cds``            1                                                    convection-diffusion solution object
``cht``            1                                                    whether the problem contains conjugate heat transfer
``dim``            1                                                    spatial dimension of ``nrs->mesh``
``dt``             3                                                    time step for previous 3 time steps
``fieldOffset``    1                                                    offset in flow solution arrays to access new component
``FU``             ``NVfields * nEXT * fieldOffset`` :math:`\checkmark` source term for each momentum equation for each step in the time stencil
``isOutputStep``   1                                                    if an output file is written on this time step
``lastStep``       1                                                    if this time step is the last time step of the run
``mesh``           1                                                    mesh used for the flow simulation
``nEXT``           1                                                    number of time steps in the time derivative stencil
``NiterU``         1                                                    number of iterations taken in last velocity solve
``NiterP``         1                                                    number of iterations taken in last pressure solve
``Nlocal``         1                                                    number of quadrature points local to this process
``Nscalar``        1                                                    number of passive scalars to solve for
``NTfields``       1                                                    number of flow-related fields to solve for (:math:`\vec{V}` plus :math:`T`)
``NVfields``       1                                                    number of velocity fields to solve for
``o_mue``          ``fieldOffset``                   :math:`\checkmark` total dynamic viscosity (laminar plus turbulent) for the momentum equation
``options``        1                                                    object containing user settings from ``.par`` file
``o_rho``          ``fieldOffset``                   :math:`\checkmark` density for the momentum equation
``P``              ``fieldOffset``                   :math:`\checkmark` pressure solution for most recent time step
``prop``           ``2 * fieldOffset``               :math:`\checkmark` total dynamic viscosity (laminar plus turbulent) and density (in this order) for the momentum equation
``U``              ``NVfields * fieldOffset``        :math:`\checkmark` velocity solution for all components for most recent time step
================== ================================= ================== ======================================================================================================

Passive Scalar Solution Fields and Simulation Settings