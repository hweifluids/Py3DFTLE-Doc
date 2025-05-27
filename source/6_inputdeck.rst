.. _inputdeck:

Input Deck for Solver Parameters
================================

To become

========================= ======================================== ============================== ========
Key                       Description                               Value(s)                        Methods
========================= ======================================== ============================== ========
``vtk_folder``             Path to folder containing VTK data files  ``<path>``                       NA
``output_dir``             Path to directory for saving FTLE results  ``<path>``                       NA
``auto_wall``              Flag for automatic wall detection (0 = off, 1 = on)  ``0`` or ``1``                 NA
``wall_xrange``            X-direction wall bounds [xmin, xmax]     ``[<float>, <float>]``         NA
``wall_yrange``            Y-direction wall bounds [ymin, ymax]     ``[<float>, <float>]``         NA
``wall_zrange``            Z-direction wall bounds [zmin, zmax]     ``[<float>, <float>]``         NA
``x``                      Number of seed points in X direction     ``<int>``                       NA
``y``                      Number of seed points in Y direction     ``<int>``                       NA
``z``                      Number of seed points in Z direction     ``<int>``                       NA
``T_range``                Integration time range [Tstart, Tend]    ``[<float>, <float>]``         NA
``dt``                     Time step for particle advection         ``<float>``                     NA
``if_backward``            Backward integration flag               ``0`` or ``1``       NA
``advc_method``            Advection scheme                        ``<string>``                   NA
``interp_method``          Interpolation method                    ``<string>``                   NA
``wall_treatment``         Wall-treatment strategy                  ``<string>``                   NA
``grad_order``             Order of spatial gradient approximation  ``<int>``                       NA
``Eigen_method``           Eigenvalue computation method           ``<string>``                   NA
``if_visual``              Visualization flag (0 = off; 1 = on)    ``0`` or ``1``                 NA
``dyn_mode``               Dynamic LCS mode selector                ``<int>``                       NA
``win_size``               Window size for dynamic LCS              ``<int> or None``              NA
``win_step``               Window step for dynamic LCS              ``<int> or None``              NA
========================= ======================================== ============================== ========

.. _flow_vars:

Flow So