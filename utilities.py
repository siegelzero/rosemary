# Miscellaneous utilities to make life easier
# Kenneth Brown 5/03/11

from functools import wraps

################################################################################
# memoization wrapper
################################################################################

def cached_function(func):
    cache = {}
    @wraps(func)
    def wrap(*args):
        if args not in cache:
            cache[args] = func(*args)
        return cache[args]
    return wrap

################################################################################

def prod(L):
    """
    prod(L):
    This returns the product of the elements of L.
    """
    p = 1
    for e in L:
        p *= e
    return p

################################################################################
