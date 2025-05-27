.. _inputdeck:

Input Deck for Solver Parameters
================================

To become

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

.. _flow_vars:

Flow So