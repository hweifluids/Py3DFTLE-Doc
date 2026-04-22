# Streamcenter+ Documentation

This repository contains the Sphinx documentation sources for `Streamcenter+`.
The documentation corresponds to the current mainline project, whose workflow is centered on a
native C++/CUDA solver and a Qt/VTK desktop GUI, while the earlier Python implementation is
retained only as a legacy archive.

## How to build locally

If you are developing documentation and want to preview the website before a pull request is
merged,
`make html` builds the user documentation as a set of interlinked HTML 
and image files.  The top-level webpage is `build/html/index.html`. To view this documentation
as a navigable HTML web page, simply navigate to the `build/html/index.html` file in
your file system and open with a web browser.
  
Note: This requires the [Sphinx](https://pypi.python.org/pypi/Sphinx) and
[sphinx_rtd_theme](https://pypi.python.org/pypi/sphinx_rtd_theme) Python packages.  Both are
available from the [Python Package Index](http://www.sphinx-doc.org://pypi.python.org/pypi).  

## How to contribute

Please create a fork of the repository and make pull/merge requests. Keep in 
mind that the number of binary files should be kept minimal. The Makefile should be 
adapted to any special build requirements.

New issues or requests are welcome to be reported.
