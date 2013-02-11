# Number Theoretic Functions
# Kenneth Brown

from bisect import bisect_left
from collections import defaultdict
from math import log, sqrt
from random import randint

from rosemary.utilities import cached_function
from rosemary.number_theory.prime_list import PRIME_LIST

import itertools

def gcd(x, y=None):
    """
    Returns the greatest common divisor of x and y.

    Returns the greatest common divisor of the integers x and y. If y is
    omitted and x is a list or tuple, then the greatest common divisor of the
    elements of x is returned.

    Input:
        * x, y - Two integers.
        * x - A list or tuple of integers.

    Output:
        * d - The integer d = gcd(x, y).

    Details:
        This function uses the classical Euclidean algorithm. If x and y are
        each of order of magnitude N, it can be shown that the overall
        complexity of the algorithm is O(log(N)^2). This is essentially the
        square of the number of digits in an operand. See section 4.2 of
        "Algorithmic Number Theory I - Efficient Algorithms"  by Bach and
        Shallit for details.

    Examples:
        >>> gcd(39064, 844)
        4
        >>> gcd([24, 18, 54])
        6
    """
    # First take care of the case where x is a list or tuple
    if isinstance(x, (list, tuple)):
        if len(x) == 1:
            return x[0]
        d = gcd(x[0], x[1])
        for i in xrange(2, len(x)):
            d = gcd(d, x[i])
        return d
    else:
        # Use the standard Euclidean algorithm
        x = abs(x)
        y = abs(y)
        # Swap x and y so that x >= y
        if x < y:
            (x, y) = (y, x)
        while y > 0:
            (x, y) = (y, x % y)
        return x

def ext_gcd(x, y):
    """
    Returns the extended gcd of x and y.

    For integers x, y with x >= y >= 0 and x > 0, this algorithm returns the
    extended gcd of x and y; i.e. the integer triple (a, b, g) such that a*x +
    b*y = g = gcd(x, y).

    Input:
        * x, y - Nonnegative integers.

    Output:
        * (a, b, g) - A tuple of integers such that a*x + b*y = g = gcd(x, y).

    Details:
        This function uses the classical Extended Euclidean algorithm. The
        complexity is essentially the same as the classical algorithm.

    Examples:
        >>> ext_gcd(5, 7)
        (3, -2, 1)
        >>> ext_gcd(25, 15)
        (-1, 2, 5)
        >>> -1 * 25 + 2 * 15
        5
    """
    if x < y:
        x, y = y, x
    (a, b, g, u, v, w) = (1, 0, x, 0, 1, y)
    while w > 0:
        q = g // w
        (a, b, g, u, v, w) = (u, v, w, a - q*u, b - q*v, g - q*w)
    return (a, b, g)

def inverse_mod(a, m):
    """
    Returns the inverse of a modulo m.

    For integers a, m with gcd(a, m) = 1, this algorithm returns the inverse of
    a modulo m; i.e. returns b such that a*b = 1 (mod m).

    Input:
        * a, m - Relatively prime integers.

    Output:
        * b - Integer such that a*b = 1 (mod m).

    Details:
        This function computes the modular inverse using the extended Euclidean
        algorithm.

    Examples:
        >>> inverse_mod(5, 17)
        7
        >>> (5 * 7) % 17
        1
    """
    e = ext_gcd(a, m)
    if e[2] != 1:
        raise ValueError('Integers not relatively prime')
    v = e[0]
    if v < 0:
        v += m
    return v

def lcm(x, y=None):
    """
    Returns the least common multiple of x and y.

    For integers x, y, this returns the least common multiple of x and y. If y
    is omitted and x is a list or a tuple, this returns the least common
    multiple of the elements on x.

    Input:
        * x, y - two integers
        * x" - list or tuple of integers

    Output:
        * m - Integer such that m = lcm(x, y).
    
    Details:
        This function uses the identity x * y  = gcd(x, y) * lcm(x, y).

    Examples:
        >>> lcm(2, 5)
        10
        >>> lcm([2, 3, 4, 5])
        60
    """
    if isinstance(x, (list, tuple)):
        if len(x) == 1:
            return x[0]
        l = lcm(x[0], x[1])
        for i in xrange(2, len(x)):
            l = lcm(l, x[i])
    else:
        l = (x * y) // gcd(x, y)
    return l

def bit_count(n):
    """
    Returns the number of set bits

    Given an integer n, this function returns the number of set bits in the
    binary expansion of n; i.e. the number of 1s appearing in the binary
    string.

    Input:
        * n - An integer.

    Output:
        * d - An integer counting the number of set bits of n.

    Details:
        This function repeatedly replaces n with n & (n - 1). Each iteration of
        this zeros out the least significant nonzero bit of n. We simply count
        the number of iterations required to annihilate n.

    Examples:
        >>> bit_count(5)
        2
        >>> bin(5)
        '0b101'
        >>> bit_count(2**8 - 1)
        8
        >>> bin(2**8 - 1)
        '0b11111111'
    """
    count = 0
    while n:
        n = n & (n - 1)
        count += 1
    return count

def integer_log(a, b):
    """
    Returns the base-b integr logarithm of a.

    Given positive integers a and b, this function returns the integer n such
    that b**n <= a < b**(n + 1).

    Input:
        * a - A positive integer. This is the argument of the logarithm.
        * b - A positive integer. This is the base of the logarithm.

    Output:
        * n - The integer with b**n <= a < b**(n + 1).

    Details:
        This function computes powers b, b**2, b**4, b**8, b**16, ..., and then
        performs a binary search. This gives the algorithm a logarithmic
        runtime instead of the linear runtime of the naive search method.

    Examples:
        >>> integer_log(8, 2)
        3
        >>> 2**3
        8
        >>> integer_log(10, 3)
        2
        >>> 3**2
        9
        >>> 3**(2 + 1)
        27
    """
    # Take care of the silly case.
    if a < b:
        return 0
    p = b
    hi = 1
    # Look at b, b^2, b^4, b^8,... to find in which interval a lives
    while p <= a:
        p = p**2
        lo = hi
        hi *= 2
    # Now we know that b^lo <= a <= b^hi perform a binary search on this
    # interval to find the exact value n so that b^n <= a < b^(n + 1)
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if b**mid > a:
            hi = mid
        else:
            lo = mid
    return lo

def integer_nth_root(n, k):
    """
    Returns the integer part of the nth root of k

    Given positive integers k and n, this returns the integer part of the nth
    root of k. This is the number r such that r**n <= k < r**(n + 1)

    Input:
        * n - A positive integer. This is the radicand.
        * k - A positive integer. This is the exponent.

    Output:
        * r - The positive integer such that r**n <= k < r**(n + 1).

    Details:
        This algorithm is based off of Newton's Method.

    Examples:
        >>> integer_nth_root(2, 10)
        3
        >>> integer_nth_root(10, 2**10)
        2
    """
    if k == 0:
        return 1
    if k == 1:
        return n

    # Here, bn is the number of bits in the base-n representation of k.
    bn = integer_log(k, n) + 1
    x = 2**(bn // 2 + 1)

    # This is just Newton's method using integer division.
    while True:
        y = ((n - 1) * x + k // (x**(n - 1))) // n
        if y >= x:
            return x
        x = y

def integer_sqrt(n):
    """
    Returns the integer part of the square root of n.,

    Given a positive integer n, this returns the integer part of the square
    root of k; i.e. the integer r with r**2 <= n < (r + 1)**2.

    Input:
        * n - A positive integer.

    Output:
        * r - The positive integer r such that r**2 <= n < (r + 1)**2.

    Details:
        The algorithm used is based on Newton's Method.

    Examples:
        >>> integer_sqrt(10)
        3
        >>> integer_sqrt(59649589127497217)
        244232653
    """
    if n in (0, 1):
        return n

    # if n is small enough, compute the square root and round
    if n <= 2**100:
        return int(sqrt(n))

    # Here, bn is the number of bits in the binary representation of n.
    bn = integer_log(n, 2) + 1
    x = 2**(bn // 2 + 1)

    # This is just Newton's method using integer division.
    while True:
        y = (x + n // x) // 2
        if y >= x:
            return x
        x = y

def is_power(x, n=None):
    """
    Determines if x in a perfect power.

    Given a positive integer x, this function returns a tuple (True, n, r) if x
    = r^n for some positive integers r, n, and returns (False, 1, x) otherwise.
    If n is not None, the function returns (True, n, r) if x = r^n for some
    positive integer r, and returns (False, 1, x) otherwise.

    Input:
        * x, n - Positive integers.
        * x - A positive integer.

    Output:
        * (b, n, r) - A triple where b is a Boolean value, and n, r are
        positive integers.

    Details:
        The algorithms is simple: if n is not None, this computes the integer
        nth root of x, and checks if the nth power of this is x. If n is none,
        the function performs trial division to find small divisors <= sqrt(x)
        and checks if x is a power of each.

    Examples:
        >>> is_power(36, 2)
        2
        >>> is_power(1331)
        3
        >>> is_power(20)
        0
    """
    if n is not None:
        r = integer_nth_root(n, x)
        if r**n == x:
            return (True, n, r)
        else:
            return (False, 1, x)
    
    d = 2
    while d*d <= x:
        m = x
        k = 0
        while m % d == 0:
            m = m // d
            k += 1
        if m == 1:
            return (True, k, d)
        d += 1
    return (False, 1, x)

def is_primitive_root(a, n):
    phi_n = euler_phi(n)
    phi_fac = factor(phi_n)

    for (p, _) in phi_fac:
        pp = power_mod(a, phi_n // p, n)
        if pp == 1:
            return False

    return True

def is_square(n):
    """
    Determines if n is a perfect square.

    Given a nonnegative integer n, this returns True of n is a perfect square,
    other returns False.

    Input:
        * n - A nonnegative integer.

    Output:
        * b - A Boolean value.

    Details:
        If n is a nonnegative integer, the function computes the integer square
        root of n, and checks if this root squared is n. If n is in factored
        form, this checks if the exponents in the prime factorization are all
        even.

    Examples:
        >>> is_square(16)
        True
        >>> is_square(15)
        False
    """
    if isinstance(n, list):
        if n[0][0] == -1:
            return False
        n_fac = n[1:]
        for (p, e) in n_fac:
            if e % 2 == 1:
                return False
        return True
    else:
        s = integer_sqrt(n)
        if s*s == n:
            return True
        return False

def is_squarefree(n):
    """
    Determines if n is squarefree.

    Given a nonnegative integer n, this return True iff n is not divisible by
    the square of an integer > 1.

    Input:
        * n - A nonnegative integer.

    Output:
        * b - A Boolean value.

    Details:
        If n is a nonnegative integer, this factors n and checks if n is
        divisible by the square of a prime.  If n is in factored form, this
        directly checks the prime factorization.

    Examples:
        >>> is_squarefree(35)
        True
        >>> is_squarefree(100)
        False
    """
    if isinstance(n, (int, long)):
        n_fac = factor(abs(n))

    elif isinstance(n, list):
        if n[0][0] == -1:
            n_fac = n[1:]
        else:
            n_fac = n
    else:
        raise ValueError("Input must be an integer or a factorization")

    for (p, e) in n_fac:
        if e > 1:
            return False

    return True

def power_mod(a, k, m):
    """
    Returns a^k (mod m).

    Given a nonnegative integers a, k and a positive integers m, this returns
    a^k (mod m).

    Input:
        * a, k, m - Integers a >= 0, b >= 0, k > 0.

    Output:
        * r - A nonnegative integer.

    Details:
        This computes a^k using a binary exponentiation method, reducing modulo
        m at each step along the way.

    Examples:
        >>> power_mod(2, 45, 17)
        15
        >>> power_mod(3, 100, 101)
        1
    """
    if a == m:
        return 0
    if k == 0:
        return 1

    t = 1
    while k:
        if k % 2:
            t = a*t % m
        k = k // 2
        a = a * a % m

    return t

def primitive_root(n):
    if n in (2, 4):
        return n - 1
    elif n % 4 == 0:
        raise ValueError("Primite Roots do not exist modulo n")

    n_fac = factor(n)

    if len(n_fac) >= 3:
        raise ValueError("Primite Roots do not exist modulo n")
    if n % 2 == 1 and len(n_fac) == 2:
        raise ValueError("Primite Roots do not exist modulo n")

    phi_n = euler_phi(n_fac)
    phi_fac = factor(phi_n)
    ll = len(phi_fac)

    for a in xrange(2, n):
        flag = True
        k = 0
        while flag and k < ll:
            p = phi_fac[k][0]
            flag = (power_mod(a, phi_n // p, n) != 1)
            k += 1
        if flag:
            return a

def carmichael_lambda(n):
    """
    carmichael_lambda(n):
    Returns the exponent of the multiplicative group of residues modulo n; i.e.
    the smallest positive integer m such that a^m = 1 (mod n) for all integers 
    a coprime to n.

    Examples:
    >>> carmichael_lambda(100)
    20
    >>> carmichael_lambda(113)
    112
    """
    if isinstance(n, (int, long)):
        n_fac = factor(abs(n))

    elif isinstance(n, list):
        if n[0][0] == -1:
            n_fac = n[1:]
        else:
            n_fac = n
    else:
        raise ValueError("Input must be an integer or a factorization")

    L = []
    for (p, e) in n_fac:
        if p != 2:
            L += [ p**(e - 1) * (p - 1) ]
        else:
            if e == 2:
                L += [ 2 ]
            elif e >= 3:
                L += [ 2**(e - 2) ]

    return lcm(L)

def euler_phi(n):
    """
    euler_phi(n):
    Returns euler's totient function of n; i.e., returns the number of positive
    integers <= n that are coprime to n.

    Examples:
    >>> euler_phi(100)
    40
    >>> euler_phi(11213)
    11212
    """
    if isinstance(n, (int, long)):
        if n == 1:
            return 1
        elif n == 0:
            raise ValueError("Zero argument in an arithmetic function")
        n_fac = factor(abs(n))

    elif isinstance(n, list):
        if n[0][0] == -1:
            n_fac = n[1:]
        else:
            n_fac = n
    else:
        raise ValueError("Input must be an integer or a factorization")

    pp = 1
    for (p, e) in n_fac:
        pp *= p**(e - 1) * (p - 1)
    return pp

def euler_phi_list(n):
    """
    Returns a list of values of euler_phi(k), for 1 <= k <= n.
    """
    L = [1] * (n + 1)
    p_list = primes(n)

    for p in p_list:
        for mul in xrange(p, n + 1, p):
            L[mul] *= (p - 1)
        pk = p*p
        while pk <= n:
            for mul in xrange(pk, n + 1, pk):
                L[mul] *= p
            pk *= p
    return L

def factorial(n):
    """
    factorial(n):
    Given a non-negative integer n, this returns n! = n*(n - 1)*(n - 2)...2*1
    The algorithm used is a binary splitting method.

    Examples:
    >>> factorial(10)
    3628800
    >>> factorial(40)
    815915283247897734345611269596115894272000000000
    """
    D = {0:1, 1:1, 2:2, 3:6, 4:24}
    if n < 0:
        raise ValueError
    elif n <= 4:
        return D[n]

    # bn is the number of binary digits of n
    pp = 1
    k = 1

    for k in xrange(1, n.bit_length()):
        lower = n >> k
        upper = n >> (k - 1)

        j = lower + 1
        if j % 2 == 0:
            j += 1

        t_pp = 1
        while j <= upper:
            t_pp *= j
            j += 2

        pp *= (t_pp**k)

    bc = bit_count(n)

    return 2**(n - bc) * pp

def moebius(n):
    """
    moebius(n):
    Given a positive integer n, this returns the value of the Moebius function
    mu of n, the multiplicative function satisfying:
        mu(1) = 1
        mu(p) = -1
        mu(p**e) = 0 for e > 1
        mu(p1*p2*...*pk) = (-1)**k

    Examples:
    >>> moebius(35)
    1
    >>> moebius(100)
    0
    >>> moebius(2)
    -1
    """
    if isinstance(n, (int, long)):
        if n == 1:
            return 1
        elif n == 0:
            raise ValueError("Zero argument in an arithmetic function")
        n_fac = factor(abs(n))

    elif isinstance(n, list):
        if n[0][0] == -1:
            n_fac = n[1:]
        else:
            n_fac = n
    else:
        raise ValueError("Input must be an integer or a factorization")

    if is_squarefree(n_fac):
        v = len(n_fac)
        return (-1)**v
    else:
        return 0

def primorial(n):
    """
    primorial(n):
    This returns the product of the first n primes p1, p2, ..., pn.
    """
    pp = 1
    N = len(PRIME_LIST)
    if n >= N:
        raise ValueError('n too large')
    for i in xrange(n):
        pp *= PRIME_LIST[i]
    return pp

def sigma(n, k = 1):
    """
    sigma(n, k = 1):
    Given integers n and k, this returns the sum of th k-th powers of the
    divisors of n.

    Examples:
    >>> sigma(9)
    13
    >>> sigma(10, 2)
    130
    """
    if isinstance(n, (int, long)):
        if n == 1:
            return 1
        elif n == 0:
            raise ValueError("Zero argument in an arithmetic function")
        n_fac = factor(abs(n))

    elif isinstance(n, list):
        if n[0][0] == -1:
            n_fac = n[1:]
        else:
            n_fac = n
    else:
        raise ValueError("Input must be an integer or a factorization")

    pp = 1
    for (p, e) in n_fac:
        pk = p**k
        pp *= (pk**(e + 1) - 1) // (pk - 1)

    return pp

def tau(n):
    """
    tau(n):
    Given an integer n, this returns the number of divisors of n.
    """
    if isinstance(n, (int, long)):
        if n == 1:
            return 1
        elif n == 0:
            raise ValueError("Zero argument in an arithmetic function")
        n_fac = factor(abs(n))

    elif isinstance(n, list):
        if n[0][0] == -1:
            n_fac = n[1:]
        else:
            n_fac = n
    else:
        raise ValueError("Input must be an integer or a factorization")

    pp = 1
    for (p, e) in n_fac:
        pp *= (e + 1)

    return pp

def factor(n):
    """
    factorization(n):
    Returns the factorization of n. Currently, this routine calls a number of
    different algorithms.

    Examples:
    >>> factor(100)
    [(2, 2), (5, 2)]
    >>> factor(537869)
    [(37, 1), (14537, 1)]
    """
    D = defaultdict(int)

    if n == 0:
        raise ValueError("Prime factorization of 0 not defined")

    # Take care of the sign
    if n < 0:
        D[-1] = 1
        n = -1 * n

    # First, strip off all small factors on n
    for p in PRIME_LIST:
        if n == 1:
            break
        elif p*p > n:
            D[n] += 1
            n = 1
            break

        while n % p == 0:
            D[p] += 1
            n = n // p

    # Next, use pollard rho
    while n > 1:
        if is_probable_prime(n):
            D[n] += 1
            n = 1
        else:
            d = factor_pollard_rho_brent(n)
            while not is_probable_prime(d):
                d = factor_pollard_rho_brent(d)

            while n % d == 0:
                D[d] += 1
                n = n // d

    p_list = sorted([ (p, D[p]) for p in D ])
    return p_list

def factor_back(F):
    """
    factor_back(F):
    Given a factorization F, this return the factored integer.

    Examples:
    >>> factor_back([(2, 2), (5, 2)])
    100
    """
    if not isinstance(F, list):
        raise ValueError("Not a factorization in factor_back")

    pp = 1
    for (p, e) in F:
        pp *= p**e

    return pp

def divisors(n):
    """
    divisors(n):
    Returns a sorted list of the positive integer divisors of n. The argument
    n can be an integer, or the factorization of an integer.

    Examples:
    >>> divisors(100)
    [1, 2, 4, 5, 10, 20, 25, 50, 100]
    >>> divisors([(2, 2), (5, 2)])
    [1, 2, 4, 5, 10, 20, 25, 50, 100]
    >>> divisors(426497)
    [1, 71, 6007, 426497]
    """
    if isinstance(n, (int, long)):
        n_fac = factor(abs(n))
    elif isinstance(n, list):
        if n[0][0] == -1:
            n_fac = n[1:]
        else:
            n_fac = n
    else:
        raise ValueError("Input must be an integer or a factorization")

    p_divs = [ p for (p, e) in n_fac ]
    div_list = []
    iter_list = ( xrange(e + 1) for (p, e) in n_fac )

    for tup in itertools.product(*iter_list):
        pp = 1
        for i in xrange(len(tup)):
            pp *= p_divs[i]**tup[i]
        div_list += [ pp ]

    div_list.sort()
    return div_list

def prime_divisors(n):
    """
    prime_divisors(n):
    Returns a list of the primes dividing n

    Examples:
    >>> prime_divisors(120)
    [2, 3, 5]
    >>> prime_divisors(5272)
    [2, 659]
    """
    if isinstance(n, (int, long)):
        n_fac = factor(abs(n))
    elif isinstance(n, list):
        if n[0][0] == -1:
            n_fac = n[1:]
        else:
            n_fac = n
    else:
        raise ValueError("Input must be an integer or a factorization")

    p_divs = [ p for (p, e) in n_fac ]
    return p_divs

def xdivisors(n):
    """
    xdivisors(n):
    Returns an iterator over the positive integer divisors of n.
    The divisors are not yielded in increasing order.
    """
    if isinstance(n, (int, long)):
        n_fac = factor(n)
    elif isinstance(n, list):
        n_fac = n
    else:
        raise ValueError("Input must be an integer or a factorization")

    p_divs = [p for (p, e) in n_fac]
    iter_list = (xrange(e + 1) for (p, e) in n_fac)
    for tup in itertools.product(*iter_list):
        pp = 1
        for i in xrange(len(tup)):
            pp *= p_divs[i]**tup[i]
        yield pp
    return

def is_probable_prime(n, a=None):
    """
    Determines if n is a strong probable prime.

    Given positive integers n and a, 1 < a < n - 1, this function returns True
    if n is a strong probable prime to the base a, otherwise returns False. If a
    is omitted, then a random value is chosen.

    Input:
        * "n, a" - Positive integers, 1 < a < n - 1.
        * "n" - A positive integer.

    Output:
        * "b" - A boolean value; True if n is a strong probable prime, False
          otherwise.

    Details:
        The algorithm used is based on a method popularized by J. Selfridge. See
        algorithm 3.5.2 of "Prime Numbers: A Computational Perspective" by
        Crandall and Pomerance for details.

    Examples:
        >>> is_probable_prime(17)

    """
    if a is None:
        a = randint(1, n - 1)

    t = n - 1
    s = 0
    while t % 2 == 0:
        s += 1
        t = t // 2

    if t % 2 == 0:
        return False
    if t == 3:
        return True

    b = power_mod(a, t, n)
    if b == 1 or b == n - 1:
        return True

    for j in xrange(1, s):
        b = (b*b) % n
        if b == n - 1:
            return True
    return False

def is_prime(n):
    """
    is_prime(n):
    Returns True if n is prime, otherwise returns False
    """
    if n < 2:
        return False
    elif n <= PRIME_LIST[-1]:
        # check to see if the number is in our list of precomputed primes
        index = bisect_left(PRIME_LIST, n)
        return PRIME_LIST[index] == n
    else:
        if not is_probable_prime(n):
            return False
        # otherwise, use trial division to ensure primality
        d = trial_division(n)
        return d == n

def is_prime_miller(n):
    """
    Determines if n is prime (pending ERH).

    Given a positive integer n, this function returns True if n is a prime;
    otherwise returns False.

    Input:
        * "n" - A positive integer.

    Output:
        * "b" - A boolean value. True if n is prime, else False.


    """
    if n < 2:
        return False

    W = min(int(2*log(n)**2), n - 1)

    for a in xrange(2, W + 1):
        if not is_probable_prime(n, a):
            return False
    return True

def next_prime(n):
    """
    next_prime(n):
    This returns the smallest prime > n.

    Examples:
    >>> next_prime(100)
    101
    >>> next_prime(10**10)
    10000000019
    """
    # Silly case
    if n == 1:
        return 2

    # No even primes > 2
    n += 1
    if n % 2 == 0:
        n += 1

    # Increment by 2 until a probable prime is found.
    while not is_prime(n):
        n += 2

    return n

def next_probable_prime(n):
    """
    next_probable_prime(n):
    This returns the smallest probable prime > n.

    Examples:
    >>> next_probable_prime(100)
    101
    >>> next_probable_prime(10**10)
    10000000019
    """
    # Silly case
    if n == 1:
        return 2

    # No even primes > 2
    n += 1
    if n % 2 == 0:
        n += 1

    # Increment by 2 until a probable prime is found.
    while not is_probable_prime(n):
        n += 2

    return n

def random_prime(a, b):
    """
    Returns a random prime in [a, b].
    """
    while True:
        n = randint(a, b)
        if is_prime(n):
            return n

def chinese(L):
    """
    chinese(L):
    Given pairwise coprime integers m_1, ..., m_k and integers x_1, ..., x_k,
    this algorithm finds an integer x such that x = x_i (mod m_i) for all i.
    The input should be of the form L = [ (x_1, m_1), ..., (x_k, m_k) ].
    """
    k = len(L)
    x, m = L[0]
    i = 1
    while i < k:
        xi, mi = L[i]
        (u, v, d) = ext_gcd(m, mi)
        if d != 1:
            raise ValueError('moduli must be coprime')
        x = u * m * xi + v * mi * x
        m *= mi
        x = x % m
        i += 1

    return x

def euler_phi_inverse(m):
    """
    euler_phi_inverse(m):
    This returns a sorted list of all integers n such that euler_phi(n) = m.
    """
    # We know euler_phi(n) is even for m >= 3
    if m % 2 == 1 and m > 1:
        return []

    mfact = factor(m)
    
    if m % 2 == 0:
        twopows = [ 2**i for i in range(mfact[0][1] + 1) ]
    else:
        twopows = [ 1 ]

    D = divisors(mfact)
    P = []

    for d in D:
        if d == 1:
            continue
        if is_prime(d + 1):
            P.append(d + 1)

    S = [(1, m)]
    for p in reversed(P):
        T = []
        for (s0, s1) in S:
            if s1 == 1:
                continue
            pk = p
            d = s1 // (p - 1)
            mmod = s1 % (p - 1)
            while mmod == 0:
                if d % 2 == 0 or d == 1:
                    T.append((pk * s0, d))
                pk *= p
                mmod = d % p
                d = d // p
        S.extend(T)

    R = set([])
    for s in S:
        if s[1] in twopows:
            j = twopows.index(s[1])
            R.add(2**(j + 1) * s[0])
            if j == 0:
                R.add(s[0])

    return sorted(R)

def jacobi_symbol(a, m):
    """
    jacobi_symbol(a, m):
    Given positive odd integer m, and integer a, this algorithm returns the
    Jacobi symbol (a|m), which for m an odd prime is also the Legendre symbol.
    """
    if m % 2 == 0:
        raise ValueError('m must be odd')
    
    a = a % m
    t = 1
    while a != 0:
        while a % 2 == 0:
            a = a // 2
            if m % 8 in (3, 5):
                t = -t
        (a, m) = (m, a)
        if a % 4 == 3 and m % 4 == 3:
            t = -t
        a = a % m

    if m == 1:
        return t

    return 0

def largest_divisor(n, b):
    """
    Returns the largest divisor d of n with d <= b.

    Given a positive integer n, this returns the largest divisor d of n such
    that d <= b.

    Input:
        * "n" - A positive integer.
        * "b" - A positive integer <= n.

    Output:
        * "d" - The largest positive integer d | n with d <= b.

    Details:
        This uses a sort of "meet in the middle" approach to find the largest
        such divisor. This method gives a nice tradeoff between time and space.

    Examples:
        >>> pp = primorial(15) * primorial(10)
        >>> pp
        3978148263608934751190154300L
        >>> tau(pp)
        1889568
        >>> b = 5413374899
        >>> time max([d for d in xdivisors(pp) if d <= b ])
        CPU times: user 7.82 s, sys: 0.01 s, total: 7.83 s
        Wall time: 7.83 s
        5413363329L
        >>> time largest_divisor(pp, b)
        CPU times: user 0.01 s, sys: 0.00 s, total: 0.01 s
        Wall time: 0.01 s
        5413363329
    """
    # If b divides n, we're done.
    target = b
    if n % target == 0:
        return target

    n_fac = factor(n)
    L = []
    for (p, e) in n_fac:
        pp = p
        for _ in xrange(e):
            L.append(pp)
            pp *= p

    k = len(L) // 2
    P1 = L[:k]
    P2 = L[k:]

    # Form all of the products in each half
    prods1 = [1]
    for p in P1:
        prods1 += [x*p for x in prods1 if x*p <= b and n % (x*p) == 0]
    prods2 = [1]
    for p in P2:
        prods2 += [x*p for x in prods2 if x*p <= b and n % (x*p) == 0]

    prods1.sort()
    prods2.sort(reverse=True)

    l1 = len(prods1)
    l2 = len(prods2)
    i = 0
    j = 0
    best = 0

    # this is a meet-in-the-middle approach
    while i < l1 and j < l2:
        current = prods1[i] * prods2[j]
        if current > target:
            j += 1
        elif current < target:
            i += 1
            if current > best:
                best = current
    return best

def maximize_divisors(B):
    """
    maximize_divisors(B):
    This returns the integer k, 1 <= k <= B, with the most divisors.
    """
    p_list = []
    pp = 1
    for prime in PRIME_LIST:
        pp *= prime
        if pp <= B:
            p_list.append(prime)
        else:
            break

    num_primes = len(p_list)
    max_tau = [1, 1]
    def backtrack(last_i, last_e, total, tau):
        if tau > max_tau[0]:
            max_tau[0] = tau
            max_tau[1] = total

        if last_i < num_primes - 1:
            next_p = p_list[last_i + 1]
            k = 1
            t_pk = next_p
            while k <= last_e:
                t_total = total * t_pk
                if t_total <= B:
                    t_tau = tau * (k + 1)
                    backtrack(last_i + 1, k, t_total, t_tau)
                else:
                    break
                k += 1
                t_pk *= next_p

    k = 1
    while 2**k <= B:
        backtrack(0, k, 2**k, k + 1)
        k += 1
    return max_tau

def sqrt_mod_p(a, p):
    """
    sqrt_mod_p(a, p):
    Given an odd prime p and an integer a with (a|p) = 1, this algorithm
    returns a solution x to x^2 = a (mod p).

    Examples:
    >>> sqrt_mod_p(10, 13)
    7
    """
    #if jacobi_symbol(a, p) != 1:
    #    raise ValueError("a is not a square modulo p")
    assert jacobi_symbol(a, p) == 1, "a must be a square modulo p"

    a = a % p

    if p % 8 in (3, 7): # p = 3, 7 (mod 8)
        x = power_mod(a, (p + 1) // 4, p)

    elif p % 8 == 5: # p = 5 (mod 8)
        x = power_mod(a, (p + 3) // 8, p)
        c = x*x % p
        if c % p  != a % p:
            x = x * power_mod(2, (p - 1) // 4, p) % p

    else: #p = 1 (mod 8)
        d = randint(2, p - 1)
        while jacobi_symbol(d, p) != -1:
            d = randint(2, p - 1)

        # Write p - 1 = 2^s * t with t odd
        s = valuation(p - 1, 2)
        t = (p - 1) >> s

        A = power_mod(a, t, p)
        D = power_mod(d, t, p)
        m = 0

        for i in xrange(s):
            # if (AD^m)^(2^(s - 1 - i)) = -1 (mod p)
            if power_mod(A * D**m, 2**(s - 1 - i), p) == p - 1:
                m += 2**i

        # Now we have AD^m = 1 (mod p)
        x = power_mod(a, (t + 1) // 2, p) * power_mod(D, m // 2, p) % p

    return x

def valuation(n, p):
    """
    valuation(n, p):
    Given integers n and p > 0, this returns the highest power of p dividing n

    Examples:
    >>> valuation(12, 2)
    2
    >>> valuation(3628800, 3)
    4
    """
    if p <= 0:
        raise ValueError('p must be positive')

    val = 0
    while n % p == 0:
        n = n // p
        val += 1
    return val

def factor_pollard_rho_brent(n):
    """
    Attempts to find a nontrivial factor of n.

    Given a composite number n, this algorithm attempts to find a nontrivial
    factor of n. The method used is Brent's improvement to the Pollard-Rho
    algorithm.
    """
    x0 = randint(0, n - 1)
    (r, q, y, g, m, c) = (1, 1, x0, 1, 3, 1)

    while True:
        x = y
        for i in xrange(r):
            y = (y*y + c) % n
        k = 0
        while True:
            ys = y
            for j in xrange(min(m, r - k)):
                y = (y*y + c) % n
                q = q * abs(x - y) % n
            g = gcd(q, n)
            k += m
            if k >= r or g > 1:
                break
        
        r *= 2
        if g > 1:
            break

    if g == n:
        while True:
            ys = (ys*ys + c) % n
            g = gcd(abs(x - ys), n)
            if g > 1:
                break
    return g

def trial_division(n, b=None):
    """
    trial_division(n, b):
    This algoritm performs trial division on n with divisors <= b. The smallest
    prime dividing n is returned if found, otherwise n is returned. We have a
    list of primes to check through first.
    """
    if b is None:
        b = integer_sqrt(n) + 1

    for d in PRIME_LIST:
        if d > b:
            # no prime divisors found <= b, so return n
            return n
        if n % d == 0:
            # return any prime divisors found
            return d

    # Next we use a segmented sieve for the rest
    for p in sieve_interval(PRIME_LIST[-1], b):
        if n % p == 0:
            return p
    return n

def lucas_lehmer(p):
    """
    lucas_lehmer(p):
    Given an odd prime p, this algorithm determines whether 2^p - 1 is prime.
    """
    v = 4
    mp = 2**p - 1
    for k in xrange(1, p - 1):
        v = (v*v - 2) % mp
    return (v == 0)

def primes(n):
    """
    Returns a list of all primes <= n.

    This program uses the sieve of Eratosthenes to generate a list of all
    primes <= n.

    Input:
        * n - A positive integer.

    Output:
        * L - a list of primes.

    Details:
        This highly optimized sieve only looks at the numbers <= n that are
        coprime to 6. The code is based on some found on stackoverflow.

    Examples:
        >>> primes(100)
        [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61,
            67, 71, 73, 79, 83, 89, 97]
        >>> len(primes(10**7))
        664579
    """
    n += 1
    offset = (n%6 > 1)

    n = { 0:n,
          1:n - 1,
          2:n + 4,
          3:n + 3,
          4:n + 2,
          5:n + 1}[n % 6]

    sieve = [True] * (n // 3)
    sieve[0] = False
    sr = integer_sqrt(n)

    for i in xrange(sr // 3 + 1):
        if sieve[i]:
            k = (3*i + 1)|1
            kk = k*k
            sieve[kk // 3::2*k] = [False]*((n // 6 - kk // 6 - 1) // k + 1)
            sieve[(kk + 4*k - 2*k*(i&1)) // 3::2*k] \
                = [False]*((n // 6 - (kk + 4*k - 2*k*(i&1)) // 6 - 1) // k + 1)

    P = [2, 3] + [(3*i + 1)|1 for i in xrange(1, n // 3 - offset) if sieve[i]]
    return P

def eratosthenes(n):
    """
    Returns a list of all primes <= n.

    This program uses the sieve of Eratosthenes to generate a list of all
    primes <= n.

    Input:
        * n - A positive integer.

    Output:
        * L - a list of primes.

    Details:
        This slightly optimized sieve only looks at the odd numbers <= n. The
        implementation is included mainly for reference.

    Examples:
        >>> eratosthenes(100)
        [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61,
            67, 71, 73, 79, 83, 89, 97]
        >>> len(eratosthenes(10**7))
        664579
    """
    m = n // 2 + 1
    b_list = [True] * m
    N = integer_sqrt(m) + 1

    for i in xrange(1, N):
        if b_list[i]:
            k = 2*i + 1
            start = (k*k - 1) // 2
            how_many = (m - start) // k + ((m - start) % k > 0)
            b_list[start::k] = [False] * how_many

    P = [2] + [ 2*i + 1 for i in xrange(1, n // 2 + n%2) if b_list[i] ]
    return P

def segmented_sieve(n):
    """
    Returns an iterator over all primes <= n.

    Given a positive integer n, this finds all primes <= n by breaking up the
    interval from 1 to n into intervals of size sqrt(d) and sieving them
    separately.

    Input:
        * n - A positive integer.

    Output:
        * X - an iterator over all primes <= n.

    Details:
        The algorithm used is a standard segmented sieve. The main improvement
        is that this algorithm uses O(sqrt(n)) space, a significant improvement
        over the O(n) space used in th standard sieve of Eratosthenes.  See
        chapter 9 of "Algorithmic Number Theory I - Efficient Algorithms" by
        Bach and Shallit for details.

    Examples:
        >>> len(segmented_sieve(10**7))
        664579
        >>> X = segmented_sieve(1000)
        >>> len([p for p in X if p % 4 == 1])
        80
    """
    # compute the primes <= sqrt(n)
    delta = integer_sqrt(n)
    initial_primes = primes(delta)
    psquares = [p*p for p in initial_primes]
    r = len(initial_primes)

    for m in xrange(0, n, delta):
        d = [1] * delta
        for i in xrange(r):
            while psquares[i] < delta:
                d[psquares[i]] = 0
                psquares[i] += initial_primes[i]
            psquares[i] -= delta
        for i in xrange(delta):
            if d[i] and 1 < m + i <= n:
                yield m + i

def sieve_interval(a, b):
    """
    Returns an iterator over all primes in the interval [a, b].

    Given positive integers a and b with a*a > b, this returns an interator
    over all primes in the interval [a, b].

    Input:
        * a - A positive integer.
        * b - A positive integer.

    Output:
        * X - An iterator over all primes in [a, b].

    Details:
        This algorithm is from Section 3.2.2 of "Prime Numbers - A
        Computational Perspective" by Crandall and Pomerance. They call it the
        "Practical Eratosthenes sieve".

    Examples:
        >>> list(sieve_interval(100, 200))
        [101, 103, 107, 109, 113, 127, 131, 137, 139, 149]
        >>> len(list(sieve_interval(10**8, 2*10**8)))
        5317482
    """
    # We need both endpoints to be even
    if a % 2 == 1:
        a -= 1
    if b % 2 == 1:
        b += 1

    bd = integer_sqrt(b)
    assert a > bd

    p_list = primes(bd)
    num_primes = len(p_list)
    max_size = min((b - a) // 2, 8*bd)
    block_size = largest_divisor(b - a, max_size)

    q_list = [ -(a + 1 + p)//2 % p for p in p_list ]

    for t in xrange(a, b, 2*block_size):
        b_list = [1] * block_size
        for k in xrange(1, num_primes):
            pk = p_list[k]
            qk = q_list[k]
            diff = block_size - qk
            how_many = diff // pk + (diff % pk > 0)
            b_list[qk:block_size:pk] = [0] * how_many
            q_list[k] = (-diff) % pk
        for j in xrange(block_size):
            if b_list[j]:
                yield t + 2*j + 1

def number_of_primes_in_interval(a, b):
    """
    Returns the number of primes p, a <= p <= b.

    Given positive integers a and b with a*a > b, this returns the number of
    primes in the interval [a, b].

    Input:
        * a - A positive integer.
        * b - A positive integer.

    Output:
        * c - The number of primes in the interval [a, b].

    Details:
        This algorithm is from Section 3.2.2 of "Prime Numbers - A
        Computational Perspective" by Crandall and Pomerance. They call it the
        "Practical Eratosthenes sieve".

    Examples:
        >>> number_of_primes_in_interval(10**8, 2*10**8)
        5317482
        >>> number_of_primes_in_interval(10**9, 2*10**9)
        47389752
    """
    # We need both endpoints to be even
    if a % 2 == 1:
        a -= 1
    if b % 2 == 1:
        b += 1

    bd = integer_sqrt(b)
    assert a > bd

    p_list = primes(bd)
    num_primes = len(p_list)
    max_size = min((b - a) // 2, 8*bd)
    block_size = largest_divisor(b - a, max_size)

    q_list = [ -(a + 1 + p)//2 % p for p in p_list ]

    count = 0
    for t in xrange(a, b, 2*block_size):
        b_list = [1] * block_size
        for k in xrange(1, num_primes):
            pk = p_list[k]
            qk = q_list[k]
            diff = block_size - qk
            how_many = diff // pk + (diff % pk > 0)
            b_list[qk:block_size:pk] = [0] * how_many
            q_list[k] = (-diff) % pk
        count += sum(b_list)
    return count
