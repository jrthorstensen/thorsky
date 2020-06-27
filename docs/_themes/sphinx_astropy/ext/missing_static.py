# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
The purpose of this extension is to give a clear warning if sphinx is expecting
a static directory to be present but it isn't.
"""

from __future__ import print_function

import os
from distutils.version import LooseVersion

from sphinx import __version__

SPHINX_LT_16 = LooseVersion(__version__) < LooseVersion('1.6')
SPHINX_LT_18 = LooseVersion(__version__) < LooseVersion('1.8')


WARNING_TEMPLATE = """
Note: The static directory '{0}' was not found. This is often because it is
      empty and you are using git. If so, you don't need it, so make
      sure it isn't included in the html_static_path setting in your conf.py
      file, otherwise Sphinx may fail the build if you are turning warnings into
      errors.
"""


def static_warning(app, config=None):

    if SPHINX_LT_16:
        info = app.info
    else:
        from sphinx.util import logging
        info = logging.getLogger(__name__).info

    for directory in app.config.html_static_path:
        if not os.path.exists(directory):
            info(WARNING_TEMPLATE.format(directory))


def setup(app):
    app.connect('build-finished', static_warning)
