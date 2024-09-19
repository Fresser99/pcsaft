from typing import NamedTuple
import numpy as np


class Param(NamedTuple):
    s: np.ndarray
    e: np.ndarray
    m: np.ndarray
    k_ij: np.ndarray
