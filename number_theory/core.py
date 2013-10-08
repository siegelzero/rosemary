# Basic number theoretic routines
# Kenneth Brown

def gcd(x, y=None):
    """
    Returns the greatest common divisor of x and y.

    Input:
        * x: int
        * y: int

    Output:
        * d: int
            The integer d = gcd(x, y).

    Details:
        This function uses the classical Euclidean algorithm. If x and y are each of order of magnitude N, it can be
        shown that the overall complexity of the algorithm is O(log(N)^2). This is essentially the square of the number
        of digits in an operand. See section 4.2 of "Algorithmic Number Theory I - Efficient Algorithms"  by Bach and
        Shallit for details.

    Examples:
        >>> gcd(39064, 844)
        4
        >>> gcd([24, 18, 54])
        6
    """
    x = abs(x)
    y = abs(y)
    if x < y:
        x, y = y, x
    while y:
        x, y = y, x % y
    return x


def ext_gcd(x, y):
    """
    Returns the extended gcd of x and y.

    For integers x, y with x >= y >= 0 and x > 0, this algorithm returns the extended gcd of x and y; i.e. the integer
    triple (a, b, g) such that a*x + b*y = g = gcd(x, y).

    Input:
        * x: int
        * y: int

    Output:
        * (a, b, g): tuple
            A tuple of integers such that a*x + b*y = g = gcd(x, y).

    Details:
        This function uses the classical Extended Euclidean algorithm. The complexity is essentially the same as the
        classical algorithm.

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


def lcm(x, y=None):
    """
    Returns the least common multiple of x and y.

    Input:
        * x: int
        * y: int

    Output:
        * m: int
            Integer such that m = lcm(x, y).
    
    Details:
        This function uses the identity x * y  = gcd(x, y) * lcm(x, y).

    Examples:
        >>> lcm(2, 5)
        10
        >>> lcm([2, 3, 4, 5])
        60
    """
    return (x * y) // gcd(x, y)


def power_mod(a, k, m):
    """
    Returns a^k (mod m).

    Given a nonnegative integers a, k and a positive integers m, this returns a^k (mod m).

    Input:
        * a: int (a >= 0)
        * k: int (k >= 0)
        * m: int (m >= 1)

    Output:
        * r: int
            r is the integer [0, m) with r = a^k (mod m)

    Details:
        This computes a^k using a binary exponentiation method, reducing modulo m at each step along the way.

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


def inverse_mod(a, m):
    """
    Returns the inverse of a modulo m.

    For integers a, m with gcd(a, m) = 1, this algorithm returns the inverse of a modulo m; i.e. returns b such that a*b
    = 1 (mod m).

    Input:
        * a: int
        * m: int

    Output:
        * b: int
            An integer such that a*b = 1 (mod m).

    Details:
        This function computes the modular inverse using the extended Euclidean algorithm.

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


def bit_count(n):
    """
    Returns the number of set bits

    Given an integer n, this function returns the number of set bits in the binary expansion of n; i.e. the number of 1s
    appearing in the binary string.

    Input:
        * n: int

    Output:
        * d: int
            The number of set bits of n.

    Details:
        This function repeatedly replaces n with n & (n - 1). Each iteration of this zeros out the least significant
        nonzero bit of n. We simply count the number of iterations required to annihilate n.

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
        * a: int (a > 0)
            This is the argument of the logarithm.

        * b: int (b > 0)
            This is the base of the logarithm.

    Output:
        * n: int
            The integer with b**n <= a < b**(n + 1).

    Details:
        This function computes powers b, b**2, b**4, b**8, b**16, ..., and then performs a binary search. This gives the
        algorithm a logarithmic runtime instead of the linear runtime of the naive search method.

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


def integer_nth_root(k, n):
    """
    Returns the integer part of the kth root of n

    Given positive integers k and n, this returns the integer part of the kth root of n. This is the number r such that
    r**k <= n < r**(k + 1)

    Input:
        * k: int (k > 0)
        * n: int (n > 0)

    Output:
        * r: int
            The positive integer such that r**k <= n < r**(k + 1).

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

    Given a positive integer n, this returns the integer part of the square root of k; i.e. the integer r with r**2 <= n
    < (r + 1)**2.

    Input:
        * n: int (n >= 0)

    Output:
        * r: int
            The positive integer r such that r**2 <= n < (r + 1)**2.

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
        return int(n**(0.5))

    # Here, bn is the number of bits in the binary representation of n.
    bn = integer_log(n, 2) + 1
    x = 2**(bn//2 + 1)

    # This is just Newton's method using integer division.
    while True:
        y = (x + n//x)//2
        if y >= x:
            return x
        x = y


def is_power(x, n=None):
    """
    Determines if x in a perfect power.

    Given a positive integer x, this function returns a tuple (True, n, r) if x = r^n for some positive integers r, n,
    and returns (False, 1, x) otherwise.  If n is not None, the function returns (True, n, r) if x = r^n for some
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


