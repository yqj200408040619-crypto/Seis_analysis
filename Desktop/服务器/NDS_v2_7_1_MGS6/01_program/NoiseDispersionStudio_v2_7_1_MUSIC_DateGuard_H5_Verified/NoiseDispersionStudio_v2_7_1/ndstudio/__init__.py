"""Reusable ambient-noise dispersion desktop and batch workflow."""
__version__='2.7.1'
from pathlib import Path
import sys
CORE=Path(__file__).resolve().parent.parent/'core'
if str(CORE) not in sys.path:sys.path.insert(0,str(CORE))
