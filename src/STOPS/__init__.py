"""Supplementary TOols for Polsalt Spectropolarimetry (STOPS)"""

__all__ = ["Split", "Join", "CrossCorrelate", "Skylines", "__version__"]

from .split import Split
from .join import Join
from .cross_correlate import CrossCorrelate
from .skylines import Skylines

__version__ = "2024.09.10"
# "Production", "Prototype", "Deprecated"
__status__ = "Development"

__description__ = "Supplementary TOols for Polsalt Spectropolarimetry (STOPS)"
__author__ = "Justin Cooper"
__email__ = "justin.jb78+Masters@gmail.com"
__url__ = "https://github.com/JustinotherGitter/STOPS"
__maintainer__ = "Justin Cooper"
__license__ = "BSD-3-Clause"
__copyright__ = "2024, Justin Cooper"
