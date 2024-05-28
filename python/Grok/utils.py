import os
import numpy as np
from contextlib import redirect_stdout, redirect_stderr, contextmanager

@contextmanager
def suppress():
    with open(os.devnull, "w") as null:
        with redirect_stdout(null):
            yield

@contextmanager
def suppress_stderr():
    with open(os.devnull, "w") as null:
        with redirect_stderr(null):
            yield            

def expand_path(path):
    return os.path.expanduser(os.path.abspath(path))

def expand_paths(input_paths):
    if isinstance(input_paths, (str, bytes)):
        input_paths = [input_paths]
    return tuple(map(expand_path, input_paths))            


def overlap(A, B):
    # get range of values in A and B that overlap    
    return (
        np.max([np.min(A), np.min(B)]), 
        np.min([np.max(A), np.max(B)])
    )


    
def to_contiguous_regions(mask):
    v = np.diff(mask.astype(int))
    indices = np.hstack([
        0, 
        np.repeat(1 + np.where(v != 0)[0], 2),
        len(mask)
    ])
            
    indices = indices.reshape((-1, 2))
    if indices.size > 0:
        offset = 0 if mask[0] else 1
        return indices[offset::2]
    return indices
        