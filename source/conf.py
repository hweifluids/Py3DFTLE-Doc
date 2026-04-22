# -*- coding: utf-8 -*-

extensions = ['sphinx.ext.mathjax', 'sphinx.ext.githubpages']
numfig = True

templates_path = []
source_suffix = '.rst'
master_doc = 'index'

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
