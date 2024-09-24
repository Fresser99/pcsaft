from typing import NamedTuple
import numpy as np


class Param(NamedTuple):
    CAS:np.ndarray
    s: np.ndarray
    e: np.ndarray
    m: np.ndarray
    k_ij: np.ndarray
    CPIG: np.ndarray
    HIGTREF: np.ndarray
    type: np.ndarray
    DVAMVDW:np.ndarray
    Tg:np.ndarray
