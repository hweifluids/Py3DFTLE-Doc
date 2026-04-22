.. _command:

Command-Line Mode (Commands)
============================

Command-line mode runs the native CUDA solver without the GUI.
It is the recommended workflow for batch jobs, remote execution, scheduler submission, and any
reproducible run driven by a checked-in ``.par`` file.
The solver reads one parameter file, loads the velocity sequence described there, computes the
requested FTLE fields, and writes VTK output to ``output_dir``.


Basic Invocation
----------------

The current native solver supports one explicit command-line option:

.. code-block:: text

   streamcenterplus_cuda --par <path-to-par-file>

In the repository workflow, the preferred entry point on Windows is the helper script:

.. code-block:: powershell

   .\CUDAsolvers\V2603\run.ps1 --par C:\path\to\case.par

On Linux, the same wrapper can be used through PowerShell:

.. code-block:: bash

   pwsh ./CUDAsolvers/V2603/run.ps1 --par /path/to/case.par

If your runtime environment is already configured, the built binary may also be launched directly.
The wrapper script remains the safer option because it resolves the newest packaged solver and adds
the required runtime-library paths automatically.


Default Behavior
----------------

If ``--par`` is omitted, the solver looks for ``default.par`` in the current working directory and
prints a warning before proceeding.
This behavior is mainly useful for quick local tests.
For production or shared workflows, an explicit ``--par`` path is strongly preferred.

The command-line solver does not currently expose ad hoc overrides for individual numerical
parameters.
In other words, all run settings are expected to come from the ``.par`` file described in
:ref:`inputdeck`.


Parameter File Advice
---------------------

The ``.par`` file is a plain text file using ``key = value`` lines, with ``#`` available for full-line
or inline comments.
For command-line use, the following practices are recommended:

- prefer absolute paths for ``vtk_folder`` and ``output_dir``;
- keep one ``.par`` file per case or per numerical setup;
- use the folder structure and qualified VTK format described in :ref:`input`;
- keep filenames stable so reruns remain reproducible.

Absolute paths are recommended because the native solver uses the paths as written.
It does not reinterpret them relative to the ``.par`` file location.


Example Run
-----------

Assuming a valid input folder and a parameter file, a typical Windows command is:

.. code-block:: powershell

   .\CUDAsolvers\V2603\run.ps1 --par C:\1_Development\Streamcenter+\CUDAsolvers\V2603\tools\case.par

During execution, the solver prints the resolved parameters, reports the number of detected VTK
frames, performs forward FTLE assembly, and optionally performs backward FTLE assembly when
``if_backward = 1``.
Progress messages are written to the terminal throughout loading, advection, FTLE assembly, and
saving.


Output Files
------------

The command-line solver writes structured VTK results into ``output_dir``.

For steady mode, namely ``dyn_mode = 0``, the output base name is
``FTLE3D_NbCU`` and the solver writes:

- ``FTLE3D_NbCU.vts``.

For dynamic mode, namely ``dyn_mode = 1``, the solver writes one ``.vts`` file per valid sliding
window:

- ``FTLE3D_NbCU_00001.vts``,
- ``FTLE3D_NbCU_00002.vts``,
- and so forth.

If at least one dynamic window is produced, the solver also writes
``FTLE3D_NbCU_dynamic.pvd`` as a collection file for time-series loading.

Within each ``.vts`` output, the point-data arrays are:

- ``FTLE_forward`` always;
- ``FTLE_backward`` only when ``if_backward = 1``.


Indexing and Frame Selection
----------------------------

The frame range is controlled by ``T_range`` in the ``.par`` file.
Its indices are applied to the input sequence after the ``.vtk`` filenames have been sorted as
described in :ref:`input`.

The current native solver:

- uses zero-based frame indexing internally;
- treats the selected range as inclusive on both ends;
- clamps the requested end index to the last available frame.

Therefore, if a folder contains :math:`N` sorted frames, the last valid explicit frame index is
:math:`N-1`.


Important Notes
---------------

- ``if_visual = 1`` is intentionally ignored in the command-line solver; interactive visualization
  belongs to the GUI workflow.
- In dynamic mode, ``win_size`` and ``win_step`` are required. Omitting either one causes the run
  to stop with an error.
- The command-line solver currently reads legacy ``.vtk`` input only. Result files are written as
  ``.vts`` and, for dynamic mode, also ``.pvd``.
- A nonzero exit code indicates a fatal error such as a missing VTK folder, absent velocity arrays,
  or an invalid parameter combination.


When to Prefer Command-Line Mode
--------------------------------

Command-line mode is preferable when:

- the run should be launched on a workstation or server without opening the GUI;
- one wants a reproducible case definition in a small text ``.par`` file;
- the same case must be repeated under several numerical settings;
- the solver is called from scripts, job schedulers, or automated regression checks.

For interactive inspection of results, however, the GUI and external post-processing tools such as
``ParaView`` remain more convenient than the terminal itself.
