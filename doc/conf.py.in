# -*- coding: utf-8 -*-
#
# TRIQS documentation build configuration file

import sys
sys.path.insert(0, "@CMAKE_CURRENT_SOURCE_DIR@/sphinxext")
sys.path.insert(0, "@CMAKE_CURRENT_SOURCE_DIR@/sphinxext/numpydoc")
sys.path.insert(0, "@CMAKE_BINARY_DIR@/python")

extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.mathjax',
              'sphinx.ext.intersphinx',
              'sphinx.ext.doctest',
              'sphinx.ext.todo',
              'sphinx.ext.viewcode',
              'sphinx.ext.autosummary',
              'sphinx.ext.githubpages',
              'sphinx_autorun',
              'matplotlib.sphinxext.plot_directive',
              'nbsphinx',
              'IPython.sphinxext.ipython_console_highlighting',
              'numpydoc']

source_suffix = '.rst'

project = '@PROJECT_NAME@'
version = '@PROJECT_VERSION@'

copyright = '2011-2020'

mathjax_path = "https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/MathJax.js?config=default"
templates_path = ['@CMAKE_CURRENT_SOURCE_DIR@/_templates']

html_theme = 'triqs'
html_theme_path = ['@CMAKE_CURRENT_SOURCE_DIR@/themes']
html_show_sphinx = False
html_context = {'header_title': 'dft tools',
                'header_subtitle': 'connecting <a class="triqs" style="font-size: 12px" href="http://triqs.github.io/triqs">TRIQS</a> to DFT packages',
                'header_links': [['Install', 'install'],
                                 ['Documentation', 'documentation'],
                                 ['Tutorials', 'tutorials'],
                                 ['Issues', 'issues'],
                                 ['About DFTTools', 'about']]}
html_static_path = ['@CMAKE_CURRENT_SOURCE_DIR@/_static']
html_sidebars = {'index': ['sideb.html', 'searchbox.html']}

htmlhelp_basename = '@PROJECT_NAME@doc'

intersphinx_mapping = {'python': ('https://docs.python.org/3.8', None), 'triqslibs': ('https://triqs.github.io/triqs/latest', None), 'triqscthyb': ('https://triqs.github.io/cthyb/latest', None)}
