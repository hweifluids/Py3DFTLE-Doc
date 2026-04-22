.. _gui:

Graphical User Interface (GUI)
==============================

The current mainline GUI of ``Streamcenter+`` is a native Qt/VTK desktop front end for the native
CUDA solver.
It is intended for users who prefer an interactive workflow for setting up ``.par`` files, launching
solver runs, monitoring progress, and previewing FTLE results without typing the full command line
manually.

The interface launches the native CUDA solver, using the same ``.par``-based workflow described in
:ref:`inputdeck` and the same structured-VTK input assumptions described in :ref:`input`.


General Workflow
----------------

In typical use, the GUI is operated in the following order:

1. Select the input ``.vtk`` folder and the output folder.
2. Fill in the solver parameters on the left-side parameter pages.
3. Save or reload a ``.par`` file if desired.
4. Start the FTLE run from the toolbar.
5. Monitor progress in the run-status strip and the live console.
6. Open the generated ``.vts`` or ``.pvd`` result for visualization.
7. Export a PNG snapshot or, for time-series results, an MP4 animation.

The GUI does not replace external post-processing entirely.
For publication-quality rendering and advanced analysis, ``ParaView`` or other VTK-capable tools
are still recommended after the initial inspection.


Top Toolbar and File Menu
-------------------------

The top toolbar combines run control and result viewing.
Its main groups are:

- the file menu, opened from the left-most menu button;
- run controls for start, pause/resume, and stop;
- a result-file row for browsing and opening ``*.vts``, ``*.vtk``, or ``*.pvd`` files;
- animation controls for time-series browsing;
- view/export controls such as parallel projection, PNG capture, and MP4 export.

The file menu currently provides these actions:

- ``Open Input Folder``;
- ``Open Output Folder``;
- ``Reload *.Par File``;
- ``Generate *.Par for FTLE``;
- ``About``.

This makes the GUI suitable both as a run launcher and as a lightweight parameter-file editor.


Left Panel: Parameter Pages
---------------------------

The left side of the window is organized as a tree and stacked pages.
In the current mainline GUI, the main pages are:

- ``File IO and Range Settings``;
- ``Numerical Methods Control``;
- ``Visualization``;
- ``Dynamic LCS``.


File IO and Range Settings
~~~~~~~~~~~~~~~~~~~~~~~~~~

This page collects the basic case definition:

- input ``.vtk`` folder;
- results folder;
- automatic or manual wall range;
- particle lattice sizes in :math:`x`, :math:`y`, and :math:`z`;
- integration range ``T_range``;
- advection step ``dt``;
- physical frame step ``physical_dt``.

Both the input folder and the output folder must already exist when a run is launched.
The current GUI checks these folders before starting the solver and warns if either one is missing.

When ``Auto wall range`` is enabled, the solver uses the input grid bounds automatically.
When it is disabled, the manual wall ranges must be filled explicitly.


Numerical Methods Control
~~~~~~~~~~~~~~~~~~~~~~~~~

This page exposes the major solver options currently supported by the native backend:

- direction: double-side or forward-only;
- wall treatment;
- advection scheme;
- interpolation method;
- gradient discretization;
- Cauchy-Green eigenvalue method.

The backend label shown on this page is informative: in the current mainline it is always the
native CUDA path.


Visualization
~~~~~~~~~~~~~

This page controls how loaded FTLE results are displayed in the embedded viewer.
The current options include:

- whether to display forward FTLE, backward FTLE, or both;
- opaque versus transparent surface rendering;
- one isovalue slider for ``FTLE_forward``;
- one isovalue slider for ``FTLE_backward``.

These sliders become meaningful only after a result file containing the corresponding FTLE arrays
has been loaded.
If one of the arrays is absent, the related control becomes unavailable automatically.


Dynamic LCS
~~~~~~~~~~~

This page switches between steady and dynamic operation:

- ``Steady`` corresponds to one still FTLE solve over the chosen range;
- ``Dynamic`` enables sliding-window processing and exposes ``win_size`` and ``win_step``.

When dynamic mode is disabled, the window-size entries are cleared and ignored.
When dynamic mode is enabled, these values should be provided consistently with the finite-time
interpretation discussed in :ref:`unsteady`.


Launching a Run
---------------

Pressing the run button causes the GUI to collect the current settings, validate the required
entries, and write an auto-generated ``.par`` file into the selected output directory.
Its filename follows the form
``GUIsolve_YYMMDDhhmmss.par``.

The GUI then launches the native solver executable with that generated parameter file.
If the solver executable cannot be found, the GUI reports where it looked and instructs the user to
build or install the solver first.

The current interaction model is centered on the parameter pages and the native solver workflow,
instead of ad hoc command-line style overrides.


Run Monitoring
--------------

During execution, two feedback channels are visible:

- a compact run-status strip under the parameter pages;
- a console panel on the lower right.

The run-status strip summarizes stages such as launching, loading, running, completion, pause, or
failure.
The console shows the forwarded solver output and GUI-side messages such as generated filenames and
log-file locations.

For each GUI-launched run, the current implementation also opens a log file in the output
directory, named after the solver process ID, for example ``12345.log``.

The pause button suspends or resumes the solver process.
The stop button terminates the running process.
These controls are intended for long FTLE jobs where an interactive interruption may be necessary.


Embedded Result Viewer
----------------------

The right side of the main window contains an embedded VTK viewer.
It can open:

- ``*.vts`` result files;
- ``*.vtk`` files for general inspection;
- ``*.pvd`` collections for time-series browsing.

Once a result is loaded, the viewer shows:

- a structured-domain outline;
- FTLE isosurfaces if ``FTLE_forward`` and/or ``FTLE_backward`` arrays are present;
- an orientation axes marker;
- a small ``Streamcenter+`` watermark;
- for ``.pvd`` time series, a time label and a frame slider.

The current visualization model is iso-surface based.
When both forward and backward FTLE fields are available, the GUI can display both simultaneously
with separate colors.
The transparency toggle is useful when overlapping surfaces should remain visible.

Parallel projection can be toggled directly from the toolbar.
The viewer uses standard interactive camera manipulation through the VTK trackball camera style.


Time-Series Playback
--------------------

When a ``.pvd`` time series is loaded, the viewer enables:

- jump to first frame;
- play or stop animation;
- jump to last frame;
- direct frame selection through the on-screen slider.

The GUI sorts and follows the time information contained in the ``.pvd`` collection file.
This makes it convenient to inspect dynamic-window outputs produced by the solver.

The MP4 export function is available only in this time-series state.
If a single still ``.vts`` file is loaded instead, the GUI informs the user that MP4 export is not
available for that case.


Saving and Reusing Parameter Files
----------------------------------

The GUI can both load an existing ``.par`` file and generate a new one from the current widget
state.
This is useful for moving between GUI-driven setup and command-line execution.

In practical terms:

- ``Reload *.Par File`` imports settings from an existing parameter file into the GUI;
- ``Generate *.Par for FTLE`` writes a user-chosen parameter file without starting a run;
- pressing ``Run`` writes a fresh auto-generated ``.par`` file and immediately launches the solver.

This workflow is especially convenient when one wants to prototype a case visually and later move
the exact setup to command-line or batch execution.


PNG and MP4 Export
------------------

Two export paths are currently integrated:

- ``PNG`` saves the current viewer image as a snapshot;
- ``MP4`` exports a time-series animation by rendering frames and calling ``ffmpeg``.

PNG export requires only a loaded scene.
MP4 export additionally requires:

- a loaded ``.pvd`` time series;
- a valid output file path;
- an ``ffmpeg`` executable available on ``PATH``.

The MP4 dialog lets the user choose:

- output filename;
- frame rate;
- frame stride;
- exported frame window.

If ``ffmpeg`` is missing, the export fails with an explicit error message.


Important Notes and Limits
--------------------------

- The GUI is a front end to the native CUDA solver.
- The GUI does not create missing input or output folders automatically; they should exist before a
  run is started.
- The GUI-side result viewer can open ``.vtk`` files, but the solver input workflow itself still
  expects a folder of qualified ``.vtk`` frames as described in :ref:`input`.
- Interactive plotting through ``if_visual`` is not the mainline command-solver workflow. Result
  inspection is expected to happen through the embedded viewer or an external post-processor.
- The embedded viewer is intentionally lightweight. For more advanced clipping, transfer functions,
  camera scripting, or batch rendering, external VTK tools remain preferable.
