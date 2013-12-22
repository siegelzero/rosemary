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
