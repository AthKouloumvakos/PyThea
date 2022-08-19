# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

import os
import sys

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import sphinx_rtd_theme

sys.path.insert(0, os.path.abspath('../../'))

# -- Project information -----------------------------------------------------

project = 'PyThea'
copyright = '2022, Athanasios Kouloumvakos'
author = 'Athanasios Kouloumvakos'
version = ''
release = ''


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.napoleon',
    'sphinx.ext.mathjax',
    'sphinx_automodapi.automodapi',
    'sphinx.ext.viewcode',
    'numpydoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.doctest',
    'sphinx.ext.inheritance_diagram',
    'jupyter_sphinx'
    # 'sphinx_gallery.gen_gallery',
]

extensions += [
    'matplotlib.sphinxext.plot_directive',
]

napoleon_google_docstring = False

# Set plotly renderer to capture _repr_html_ for sphinx-gallery
try:
    import plotly.io.renderers
except ImportError:
    pass
else:
    plotly.io.renderers.default = 'sphinx_gallery'


# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_logo = '../logo/pythea_logo_wb.png'
html_theme_options = {
    'logo_only': False,
    'display_version': True,
}
# html_favicon = "./_static/icon.ico"

nitpicky = True
nitpick_ignore = [('py:class', "Unit('deg')"),
                  ('py:class', "Unit('m')"),
                  ('py:class', "Unit('rad')")]

intersphinx_mapping = {'python': ('https://docs.python.org/3/', None),
                       'xarray': ('https://xarray.pydata.org/en/stable/', None),
                       'astropy': ('https://docs.astropy.org/en/stable/', None),
                       'scipy': ('https://docs.scipy.org/doc/scipy', None),
                       'sunpy': ('https://docs.sunpy.org/en/stable/', None),
                       'matplotlib': ('https://matplotlib.org/stable/', None),
                       'numpy': ('https://numpy.org/doc/stable/', None),
                       'pandas': ('http://pandas.pydata.org/pandas-docs/dev', None)
                       }

# default_role = 'any'
# automodapi_inheritance_diagram = False
# automodsumm_inherited_members = True
