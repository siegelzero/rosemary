# Special Functions

import math

prec = 1e-14

################################################################################
# Special constants
################################################################################

euler_gamma = 0.57721566490153286060651209008240243104
pi = 3.1415926535897932384626433832795028842

################################################################################

def logarithmic_integral(x):
    """
    logarithmic_integral(x):
    This returns the logarithmic integral of x; i.e.
    \int_{0}^{x} 1 / \log(t) \, dt
    """
    log_x = math.log(x)
    s = euler_gamma + math.log(log_x)

    k = 1
    num = log_x
    den = k
    while True:
        t = num / (den * k)
        s += t
        if abs(t) < prec:
            return s
        k += 1
        num *= log_x
        den *= k

################################################################################
