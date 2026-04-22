# -*- coding: utf-8 -*-

extensions = ['sphinx.ext.mathjax', 'sphinx.ext.githubpages']
numfig = True

templates_path = []
source_suffix = '.rst'
master_doc = 'index'

# General information about the project.
project = u'Streamcenter+'
copyright = u'2025, Huanxia Wei'
author = u'Huanxia Wei and his buddies'

version = u'0.1'
release = u'0.1'

language = 'en'
exclude_patterns = []
pygments_style = 'sphinx'
todo_include_todos = False

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

htmlhelp_basename = 'StreamcenterPlus-Doc'

latex_elements = {}

# Output file base name for HTML help builder.
htmlhelp_basename = 'StreamcenterPlus-Doc'


# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',

    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
#latex_documents = [
#    (master_doc, 'nekRS.tex', u'nekRS Documentation',
#     u'Paul Fischer, \\and James Lottes, \\and Stefan Kerkemeier, \\and Oana Marin, \\and Katherine Heisey, \\and Aleks Obabko, \\and Elia Merzari, \\and Yulia Peet, \\and Ronald Rahaman', 'manual'),
#]


# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'StreamcenterPlus', u'Streamcenter+ Documentation',
     [author], 1)
]

texinfo_documents = [
    (master_doc, 'StreamcenterPlus', u'Streamcenter+ Documentation',
     author, 'StreamcenterPlus', 'A native CUDA and Qt/VTK toolbox for FTLE and LCS research.',
     'Miscellaneous'),
]

epub_title = project
epub_author = author
epub_publisher = author
epub_copyright = copyright

epub_exclude_files = ['search.html']
