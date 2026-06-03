# =============================================================================
# load_libraries_python.py — Python dependencies for pseudotime analysis
# =============================================================================
# Loaded automatically by capitulo3_pseudotime.py.
# Do not modify unless adding new packages.

import os
import sys
import warnings

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from sklearn.preprocessing import scale

import scanpy as sc
import scFates as scf
import palantir

import matplotlib

def _running_in_jupyter():
    try:
        from IPython import get_ipython
        shell = get_ipython()
        return shell is not None and "IPKernelApp" in shell.config
    except Exception:
        return False

if _running_in_jupyter():
    try:
        from IPython import get_ipython
        get_ipython().run_line_magic("matplotlib", "inline")
    except Exception:
        pass
else:
    matplotlib.use("Agg")  # non-interactive backend for script execution

import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings("ignore")
sc.settings.verbosity = 3
sc.settings.logfile = sys.stdout
