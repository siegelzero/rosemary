from rosemary.number_theory.elementary import integer_sqrt, inverse_mod, power_mod

def discrete_log_bf(a, b, n):
    """
    finds the least positive integer k such that a**k = b (mod p) by direct
    check
    """
    if b == 1:
        return 0

    ak = a
    for k in xrange(1, n):
        if ak == b:
            return k
        ak = (ak * a) % n

    # We've tried all possibilities at this point.
    raise ValueError("No Solution")


def discrete_log_shanks(a, b, n):
    """
    finds the least positive integer k such that a**k = b (mod p) by direct
    check
    """
    m = integer_sqrt(n)
    ak = 1
    D = {0:1}
    for k in xrange(1, m):
        ak = (ak * a) % n
        D[ak] = k

    am = inverse_mod(ak * a, n)
    g = b
    for i in xrange(m):
        if g in D:
            return i * m + D[g]
        g = (g * am) % n

    raise ValueError("No Solution")

