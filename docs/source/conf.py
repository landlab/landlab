import os
import pathlib
import re
import sys
from datetime import date

import packaging
import tomli


def get_version_from_file(path):
    with open(path) as fp:
        match = re.search(r'__version__\s*=\s*[\'"]([^\'"]+)[\'"]', fp.read())
        if match:
            version = match.group(1)
        else:
            raise ValueError(f"version string not found ({path})")
    return packaging.version.Version(version)


src_dir = os.path.abspath(
    os.path.join(os.path.dirname(__file__), os.pardir, os.pardir, "src")
)
sys.path.insert(0, src_dir)
docs_dir = pathlib.Path(__file__).parent
version_file = os.path.join(src_dir, "landlab", "_version.py")

# -- General configuration -----------------------------------------------------

extensions = [
    "myst_parser",
    "nbsphinx",
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx.ext.todo",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
    "sphinx_copybutton",
    "sphinx_inline_tabs",
    "sphinxcontrib.towncrier",
    "sphinx_jinja",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# The suffix of source filenames.
source_suffix = ".rst"

# The encoding of source files.
# source_encoding = 'utf-8-sig'

# Regex for links that we know work in browser, but do not work in sphinx/CI
# (BE VERY CAREFUL ADDING LINKS TO THIS LIST)
if os.getenv("GITHUB_ACTIONS"):
    linkcheck_ignore = [
        # Added by KRB Dec 2019, at this point two links match this pattern
        r"https://pubs.geoscienceworld.org/gsa/geology.*",
        r"https://doi.org/10.1130/*",  # Added by KRB Jan 2019. Four links match this pattern
        r"https://dx.doi.org/10.1029/2011jf002181",  # Added by EWHH April 2020
        r"https://doi.org/10.1029/2019JB018596",  # Added by EWHH April 2020
        r"https://doi.org/10.3133/pp294B",  # Added by EWHH September 2021
        #     r"https://yaml.org/start.html",  # Added by EWHH September 2021
    ]
    linkcheck_retries = 5

master_doc = "index"

project = "landlab"
copyright = str(date.today().year) + ", The Landlab Team"

v = get_version_from_file(version_file)
version = f"{v.major}.{v.minor}"
release = v.public

language = "en"

# There are two options for replacing |today|: either, you set today to some
# non-false value, then it is used:
# today = ''
# Else, today_fmt is used as the format for a strftime call.
# today_fmt = '%B %d, %Y'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ["_build", "**/*.so", "**/*.pyd", "**/*.pyx", "**/*.c", "**/*.cpp"]

# The reST default role (used for this markup: `text`) to use for all documents.
# default_role = None

# If true, '()' will be appended to :func: etc. cross-reference text.
add_function_parentheses = False

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
add_module_names = False

# If true, sectionauthor and moduleauthor directives will be shown in the
# output. They are ignored by default.
show_authors = True

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"
pygments_dark_style = "monokai"

# A list of ignored prefixes for module index sorting.
modindex_common_prefix = ["landlab."]

# If true, keep warnings as "system message" paragraphs in the built documents.
# keep_warnings = False

# The default highlight language
highlight_language = "none"

# -- Options for HTML output ---------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
# html_theme = "alabaster"
html_theme = "furo"
html_title = "landlab"
html_logo = "_static/landlab_logo.png"

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
html_theme_options = {
    "announcement": "<em>Landlab 2.9 released!</em>",
    "source_repository": "https://github.com/landlab/landlab/",
    "source_branch": "master",
    "source_directory": "docs/source",
    "sidebar_hide_name": True,
    "footer_icons": [
        {
            "name": "power",
            "url": "https://csdms.colorado.edu",
            "html": """
                <svg
                  stroke="currentColor"
                  fill="currentColor"
                  stroke-width="0"
                  version="1.1"
                  viewBox="0 0 16 16"
                  height="1em"
                  width="1em"
                  xmlns="http://www.w3.org/2000/svg"
                >
                  <path
                    d="M6 0l-6 8h6l-4 8 14-10h-8l6-6z"
                  ></path>
                </svg>
                <b><i>Powered by CSDMS</i></b>
            """,
            "class": "",
        },
    ],
}

# Add any paths that contain custom themes here, relative to this directory.
# html_theme_path = []

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
# html_title = None

# A shorter title for the navigation bar.  Default is the same as html_title.
# html_short_title = None

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
# html_favicon = None
html_favicon = "_static/favicon.ico"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
# html_last_updated_fmt = '%b %d, %Y'

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
# html_use_smartypants = True

# Custom sidebar templates, maps document names to template names.
# html_sidebars = {}


# Additional templates that should be rendered to pages, maps page names to
# template names.
# html_additional_pages = {}

# If false, no module index is generated.
# html_domain_indices = True

# If false, no index is generated.
# html_use_index = True

# If true, the index is split into individual pages for each letter.
# html_split_index = False

# If true, links to the reST sources are added to the pages.
# html_show_sourcelink = True

# If true, "Created using Sphinx" is shown in the HTML footer. Default is True.
# html_show_sphinx = True

# If true, "(C) Copyright ..." is shown in the HTML footer. Default is True.
# html_show_copyright = True

# If true, an OpenSearch description file will be output, and all pages will
# contain a <link> tag referring to it.  The value of this option must be the
# base URL from which the finished HTML is served.
# html_use_opensearch = ''

# This is the file name suffix for HTML files (e.g. ".xhtml").
# html_file_suffix = None

# Output file base name for HTML help builder.
htmlhelp_basename = "landlabdoc"


todo_include_todos = True
# latex_elements = dict(preamble='\\usepackage{amsmath}')

# -- Options for napoleon extension --------------------------------------------

napoleon_numpy_docstring = True
napoleon_google_docstring = False
napoleon_include_init_with_doc = True
napoleon_include_special_with_doc = True

# -- Options for towncrier_draft extension --------------------------------------------

# or: 'sphinx-release', 'sphinx-version'
towncrier_draft_autoversion_mode = "sphinx-release"
towncrier_draft_include_empty = True
towncrier_draft_working_directory = pathlib.Path(docs_dir).parent.parent

# -- Options for intersphinx extension ---------------------------------------

# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "xarray": ("https://docs.xarray.dev/en/stable/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
}

with open("../index.toml", "rb") as fp:
    cats = tomli.load(fp)
cats["grids"].pop("ModelGrid")

jinja_contexts = {"llcats": cats}

# -- Options for autodoc

# selects what content will be inserted into the main body of an autoclass
# directive: 'class'(default), 'both', or 'init'
autoclass_content = "both"
autodoc_typehints = "description"
autodoc_class_signature = "separated"

with open(os.path.join(src_dir, "../cython-files.txt")) as fp:
    cython_files = {fname.strip() for fname in fp.readlines()}

autodoc_mock_imports = [
    "richdem",
    "bmipy",
    "importlib-resources",
    "matplotlib",
    "netcdf4",
    "pandas",
    "pyshp",
    "pyyaml",
    "rich-click",
    "scipy",
    "statsmodels",
    "xarray",
    "shapefile",
    "mpl_toolkits",
    "rich_click",
    "pylab",
] + [path[4:-4].replace("/", ".") for path in cython_files]

autodoc_default_options = {
    "maxdepth": 2,
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
    "inherited-members": "int",
    "ignore-module-all": True,  # Ignore __all__ directives in modules
}

nbsphinx_execute = "never"
nbsphinx_thumbnails = (
    {}
    if os.getenv("GITHUB_ACTIONS")
    else {
        "teaching/**/*": "_static/favicon.ico",
        "teaching/*": "_static/favicon.ico",
        "tutorials/**/*": "_static/favicon.ico",
        "tutorials/*": "_static/favicon.ico",
    }
)

# This is processed by Jinja2 and inserted before each notebook
nbsphinx_prolog = """
{% set docname = 'notebooks/' + env.doc2path(env.docname, base=None) %}

.. note::

    This page was generated from a jupyter notebook_.

.. _notebook: https://github.com/landlab/landlab/blob/{{ env.config.release|e }}/{{ docname|e }}
"""

nbsphinx_epilog = """
----

Generated by nbsphinx_ from a Jupyter_ notebook.

.. _nbsphinx: https://nbsphinx.readthedocs.io/
.. _Jupyter: https://jupyter.org/
"""

myst_enable_extensions = ["colon_fence", "deflist"]
