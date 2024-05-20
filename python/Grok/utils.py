import os
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