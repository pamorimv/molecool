"""
molecool
A Python package for analyzing and visualizing pdb and xyz files. For MolSSI May webinar series.
"""

# Add imports here
from .functions import canvas, zen
from .molecule import build_bond_list
from .visualize import draw_molecule, bond_histogram
from .measure import calculate_distance, calculate_angle
from .atom_data import atom_colors, atomic_weights
from .io import open_pdb, open_xyz, write_xyz


# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
