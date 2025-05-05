import numpy as np
import time

def generate_basis(N):
    """Generates all 2^N basis states as integers."""
    return np.arange(2**N)

def bit_parity(x):
    """Computes the parity (number of 1s mod 2) of a binary integer."""
    return bin(x).count("1") % 2

def timing_decorator(func):
    """Decorator to measure execution time of functions."""
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        print(f"{func.__name__} took {end_time - start_time:.4f} seconds")
        return result
    return wrapper