import numpy as np

def build_library():
    raise NotImplementedError

def build_optimizer():
    raise NotImplementedError

def fit_sindy(t: np.ndarray, x: np.ndarray, dt, lib, opt):
    raise NotImplementedError
