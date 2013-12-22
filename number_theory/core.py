################################################################################
# Basic number theoretic routines
################################################################################

def gcd(x, y):
    """
    Returns the greatest common divisor of x and y.

    Input:
        * x: int
        * y: int

    Output:
        * d: int
            The integer d = gcd(x, y).

    Examples:
        >>> gcd(39064, 844)
        4
        >>> gcd(7, 5)
        1

    Details:
        This function uses the classical Euclidean algorithm. If x and y are
        each of order of magnitude N, it can be shown that the overall
        complexity of the algorithm is O(log(N)^2). This is essentially the
        square of the number of digits in an operand. See section 4.2 of
        "Algorithmic Number Theory I - Efficient Algorithms"  by Bach and
        Shallit for details.
    """
    x = abs(x)
    y = abs(y)
    if x < y:
        x, y = y, x

    while y:
        x, y = y, x % y
    return x


def gcd_list(L):
    """
    Returns the greatest common divisor of the elements of L.

    Input:
        * L: list
            A list of integers.

    Output:
        * d: int
            The integer d = gcd(L[0], L[1], ..., L[n - 1]), where len(L) = n.

    Examples:
        >>> gcd_list([2, 4, 6])
        2

    Details:
        This function recursively uses the fact that
        gcd(a, b, c) = gcd(a, gcd(b, c)).
    """
    num_values = len(L)
    if num_values == 1:
        return L[0]

    partial = gcd(L[0], L[1])
    for i in xrange(2, num_values):
        partial = gcd(partial, L[i])
    return partial


def ext_gcd(x, y):
    """
    Returns the extended gcd of x and y.

    For integers x, y, this algorithm returns the extended gcd of x and y; i.e.
    the integer triple (a, b, d) such that a*x + b*y = d = gcd(x, y).

    Input:
        * x: int
        * y: int

    Output:
        * (a, b, d): tuple
            A tuple of integers such that a*x + b*y = d = gcd(x, y).

    Examples:
        >>> ext_gcd(5, 7)
        (3, -2, 1)
        >>> 3*5 - 2*7
        1
        >>> ext_gcd(12, 10)
        (1, -1, 2)
        >>> ext_gcd(12, -18)
        (-1, -1, 6)
        >>> ext_gcd(-3, -1)
        (0, -1, 1)

    Details:
        This function uses the classical Extended Euclidean algorithm. The
        complexity is essentially the same as the classical algorithm. See
        Algorithm 2.1.4 in "Prime Numbers - A Computational Perspective" by
        Crandall and Pomerance for details.
    """
    switch = False
    if x < y:
        x, y = y, x
        switch = True

    sign_x = 1 if x >= 0 else -1
    sign_y = 1 if y >= 0 else -1
    x = abs(x)
    y = abs(y)

    (a, b, d, u, v, w) = (1, 0, x, 0, 1, y)
    while w:
        q = d//w
        (a, b, d, u, v, w) = (u, v, w, a - q*u, b - q*v, d - q*w)

    if switch:
        a, b = b, a
        sign_x, sign_y = sign_y, sign_x
    return (sign_x*a, sign_y*b, d)


def lcm(x, y):
    """
    Returns the least common multiple of x and y.

    Input:
        * x: int
        * y: int

    Output:
        * m: int
            Integer such that m = lcm(x, y).

    Examples:
        >>> lcm(2, 5)
        10

    Details:
        This function uses the identity x*y  = gcd(x, y)*lcm(x, y).
    """
    x = abs(x)
    y = abs(y)
    m = (x*y)//gcd(x, y)
    return m


def lcm_list(L):
    """
    Returns the least common multiple of the elements of L.

    Input:
        * L: list
            A list of integers.

    Output:
        * l: int
            The integer l = gcd(L[0], L[1], ..., L[n - 1]), where len(L) = n.

    Examples:
        >>> lcm_list([2, 3, 4])
        12

    Details:
        This function recursively uses the fact that
        lcm(a, b, c) = lcm(a, lcm(b, c)).
    """
    num_values = len(L)
    if num_values == 1:
        return L[0]

    partial = lcm(L[0], L[1])
    for i in xrange(2, num_values):
        partial = lcm(partial, L[i])
    return partial


def power_mod(a, k, m):
    """
    Returns a^k (mod m).

    Given a nonnegative integers a, k and a positive integers m, this returns
    a^k (mod m).

    Input:
        * a: int (a >= 0)
        * k: int (k >= 0)
        * m: int (m >= 1)

    Output:
        * r: int
            r is the integer in [0, m) with r = a^k (mod m)

    Examples:
        >>> power_mod(2, 45, 17)
        15
        >>> power_mod(3, 100, 101)
        1

    Details:
        This computes a^k using a binary exponentiation method, reducing modulo
        m at each step along the way. Note that the builtin python function pow
        behaves exactly the same as this, and is typically faster, so its usage
        is preferred. See Algorithm 1.2.1 in "A Course in Computational
        Algebraic Number Theory" by Cohen for details.
    """
    if a < 0:
        raise ValueError("power_mod: Must have a >= 0.")
    if k < 0:
        raise ValueError("power_mod: Must have k >= 0.")
    if m < 1:
        raise ValueError("power_mod: Must have m >= 1.")

    if a % m == 0:
        return 0
    if k == 0:
        return 1
    if a == 0:
        return 0

    a = a % m
    t = 1
    while k:
        if k % 2:
            t = a*t % m
        k = k//2
        a = a*a % m
    return t


def inverse_mod(a, m):
    """
    Returns the inverse of a modulo m.

    For integers a, m with gcd(a, m) = 1, this algorithm returns the inverse of
    a modulo m; i.e. returns b such that a*b = 1 (mod m).

    Input:
        * a: int (a >= 1)
        * m: int (m >= 2)

    Output:
        * b: int
            An integer such that a*b = 1 (mod m).

    Examples:
        >>> inverse_mod(5, 17)
        7
        >>> 5*7 % 17
        1

    Details:
        This function computes the modular inverse using the extended Euclidean
        algorithm.
    """
    if m < 2:
        raise ValueError("inverse_mod: Must have m >= 2.")

    (x, y, d) = ext_gcd(a, m)
    if d != 1:
        raise ValueError("inverse_mod: Integers must be relatively prime.")
    return x % m


def bit_count(n):
    """
    Returns the number of set bits.

    Given an integer n, this function returns the number of set bits in the
    binary expansion of n; i.e. the number of 1s appearing in the binary
    representation of n.

    Input:
        * n: int

    Output:
        * d: int
            The number of set bits of n.

    Examples:
        >>> bit_count(5)
        2
        >>> bin(5)
        '0b101'
        >>> bit_count(2**8 - 1)
        8
        >>> bin(2**8 - 1)
        '0b11111111'

    Details:
        This function repeatedly replaces n with n & (n - 1). Each iteration of
        this zeros out the least significant nonzero bit of n. We simply count
        the number of iterations required to annihilate n.
    """
    n = abs(n)
    count = 0
    while n:
        n = n & (n - 1)
        count += 1
    return count


def integer_log(b, n):
    """
    Returns the base-b integer logarithm of n.

    Given a positive integers n and base b, this function returns the integer k
    such that b**k <= a < b**(k + 1).

    Input:
        * b: int (b >= 2)
            The base of the logarithm.

        * n: int (n >= 1)
            The argument of the logarithm.

    Output:
        * k: int
            The integer such that b**k <= a <= b**(k + 1).

    Examples:
        >>> integer_log(2, 100)
        6
        >>> integer_log(5, 30)
        2

    Details:
        This function computes powers b, b**2, b**4, b**8, b**16, ..., and then
        performs a binary search. This gives the algorithm a logarithmic runtime
        instead of the linear runtime of the naive search method.
    """
    if b <= 1:
        raise ValueError("integer_log: Must have b >= 2.")
    if n <= 0:
        raise ValueError("integer_log: Must have n >= 1.")

    if n < b:
        return 0

    p = b
    hi = 1
    # Look at b, b^2, b^4, b^8,... to find in which interval a lives
    while p <= n:
        p = p**2
        lo = hi
        hi *= 2
    # Now we know that b^lo <= a <= b^hi perform a binary search on this
    # interval to find the exact value n so that b^n <= a < b^(n + 1)
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if b**mid > n:
            hi = mid
        else:
            lo = mid
    return lo


def integer_nth_root(k, n):
    """
    Returns the integer part of the kth root of n

    Given positive integers k and n, this returns the integer part of the kth
    root of n. This is the number r such that r**k <= n < r**(k + 1)

    Input:
        * k: int (k >= 1)
        * n: int (n >= 0)

    Output:
        * r: int
            The positive integer such that r**k <= n < r**(k + 1).

    Examples:
        >>> integer_nth_root(2, 10)
        3
        >>> integer_nth_root(10, 2**10)
        2

    Details:
        This algorithm is based off of Newton's Method.
    """
    if n < 0:
        raise ValueError("integer_nth_root: Must have n >= 0")
    if k < 1:
        raise ValueError("integer_nth_root: Must have k >= 1")

    if k == 0:
        return 1
    if k == 1:
        return n
    if n in (0, 1):
        return n

    bn = integer_log(k, n) + 1
    x = 2**(bn//2 + 1)

    while True:
        y = ((k - 1)*x + n//(x**(k - 1)))//k
        if y >= x:
            return x
        x = y


def integer_sqrt(n):
    """
    Returns the integer part of the square root of n.,

    Given a positive integer n, this returns the integer part of the square root
    of k; i.e. the integer r with r**2 <= n < (r + 1)**2.

    Input:
        * n: int (n >= 0)

    Output:
        * r: int
            The positive integer r such that r**2 <= n < (r + 1)**2.

    Examples:
        >>> integer_sqrt(10)
        3
        >>> integer_sqrt(59649589127497217)
        244232653

    Details:
        The algorithm used is based on Newton's Method. See Algorithm 1.7.1 in
        'A Course in Computational Algebraic Number Theory' by Cohen for
        details.
    """
    if n < 0:
        raise ValueError("integer_sqrt: Must have n >= 0.")

    if n in (0, 1):
        return n

    # if n is small enough, compute the square root and round
    if n <= 2**60:
        return int(n**(0.5))

    # Here, bn is the number of bits in the binary representation of n.
    #bn = integer_log(n, 2) + 1
    bn = n.bit_length()
    x = 2**(bn//2 + 1)

    # This is just Newton's method using integer division.
    while True:
        y = (x + n//x)//2
        if y >= x:
            return x
        x = y


def is_power(n, k=None):
    """
    Determines if n is a perfect power.

    If n is a perfect power, this returns (b, k) where n = b^k with k maximal.
    If the optional parameter k is given, this returns (b, k) if n = b^k for
    some b. Otherwise, returns False.

    Input:
        * n: int or list (n >= 2)
            The value of n can be an int or a factorization.

        * k: int (k >= 1) (default=None)

    Output cases:
        * (b, k): tuple of ints
            Values such that n = b^k if such values exist.
        
        * Returns False if no such values exist.

    Examples:
        >>> is_power(36)
        (6, 2)
        >>> is_power(81)
        (3, 4)
        >>> is_power(81, 2)
        (9, 2)
        >>> is_power(1330)
        False
        >>> is_power([(2, 4), (5, 4)])
        (10, 4)
        >>> is_power([(2, 4), (5, 4)], 2)
        (100, 2)
        >>> is_power([(2, 3), (5, 2)])
        False

    Details:
        If n is given as an int and a value of k is given, then this function
        computes the integer kth root r of n and checks if r^k = n. In the case
        where no such k is given, this looks at each k >= 2, and performs a
        binary search for values b such that n = b^k. Details for this case can
        be found in "Primality Testing in Polynomial Time" by Dietzfelbinger. In
        particular, see Algorithm 2.3.5.

        If n is given as a factorization and a value of k is given, this
        function checks if each exponent in the prime factorization of n is
        divisible by k. If no such k is given, this computes the gcd d of the
        exponents in the prime factorization and checks if d > 1.
    """
    if k is not None and k < 1:
        raise ValueError("is_power: Must have k >= 1.")

    if isinstance(n, list):
        if n[0][0] == -1:
            raise ValueError("is_power: Must have n >= 2.")
        n_factorization = n[:]

        kth_root = 1
        if k is not None:
            for (p, e) in n_factorization:
                if e % k == 0:
                    kth_root *= p**(e//k)
                else:
                    return False
            return (kth_root, k)
        else:
            exponents = [e for (p, e) in n_factorization]
            d = gcd_list(exponents)
            if d == 1:
                return False
            else:
                kth_root = 1
                for (p, e) in n_factorization:
                    kth_root *= p**(e//d)
                return (kth_root, d)
    else:
        if n < 2:
            raise ValueError("is_power: Must have n >= 2.")

        if k is not None:
            if k == 1:
                return (n, 1)
            root = integer_nth_root(k, n)
            if root**k == n:
                return (root, k)
            return False

        k = 2
        while 2**k <= n:
            lo = 1
            hi = n
            while hi - lo > 1:
                mid = (lo + hi)//2
                power = mid**k
                if power == n:
                    # Here, we've found that n is a perfect power. However, at this
                    # point, b is the smallest exponent. We want the largest.
                    exponent = k
                    base = mid
                    while True:
                        next_step = is_power(base)
                        if next_step:
                            base = next_step[0]
                            exponent *= next_step[1]
                        else:
                            return (base, exponent)
                elif power < n:
                    lo = mid
                else:
                    hi = mid
            k += 1
        return False


def is_square(n):
    """
    Determines if n is a perfect square.

    Given a nonnegative integer n, this returns False if n is not a perfect
    square, and returns the square root of n otherwise.

    Input:
        * n: int or list (n >= 0)
            The value of n can be an int or a factorization.

    Output Cases:
        * sqrt: int
            If n is a perfect square, the square root of n is returned.

        * False if n is not a perfect square.

    Examples:
        >>> is_square(16)
        4
        >>> is_square(15)
        False
        >>> is_square([(2, 2), (5, 2)])
        10
        >>> is_square([(2, 3), (5, 2)])
        False

    Details:
        If n is a nonnegative integer, the function computes the integer square
        root of n, and checks if this root squared is n. If n is in factored
        form, this checks if the exponents in the prime factorization are all
        even.
    """
    if isinstance(n, list):
        if n[0][0] == -1:
            raise ValueError("is_square: Must have n >= 0.")
        n_factorization = n[:]
        sqrt = 1
        for (p, e) in n_factorization:
            if e % 2 == 1:
                return False
            else:
                sqrt *= p**(e//2)
        return sqrt
    else:
        if n < 0:
            raise ValueError("is_square: Must have n >= 0.")
        sqrt = integer_sqrt(n)
        if sqrt*sqrt == n:
            return sqrt
        return False


def chinese(L):
    """
    Returns a solution to the congruences given in L.

    Given a list of congruences as (residue, modulus) pairs, where all moduli
    are pairwise coprime, this returns the unique solution n in [0, M) to the
    congruences, where M is the product of the moduli.

    Input:
        * L: list
            L is a list of congruences of the form (residue, modulus).
    
    Output:
        * x: int
            x is the least positive integer satisfying all of the congruences.

    Examples:
        >>> chinese([(0, 4), (1, 25)])
        76

    Details:
        This algorithm is based on the standard constructive proof of the
        Chinese Remainder Theorem. For a single CRT application, this is the
        function to use. If solving multiple systems of congruences with the
        same set of moduli, the function 'chinese_preconditioned' should be
        used.
    """
    (x, m) = L[0]
    rest = L[1:]
    for (xi, mi) in rest:
        (u, v, d) = ext_gcd(m, mi)
        if d != 1:
            raise ValueError('chinese: Moduli must be coprime.')
        x = u*m*xi + v*mi*x
        m *= mi
        x = x % m
    return x


def crt_preconditioning_data(moduli):
    """
    Returns the preconditioning data for the set of moduli.

    Given a set of moduli, this computes the preconditioning data necessary for
    the chinese_preconditioned function.

    Input:
        * moduli: list
            A list of pairwise coprime positive integers.

    Output:
        * data: tuple
            This is a tuple containing:
                * r: The number of moduli.

                * partial_products: A list of partial products of the moduli.

                * inverse_list: A list of the inverse of each partial product,
                    modulo the next modulus in the list.

                * product: The product of the moduli.
    """
    r = len(moduli)
    partial_products = [1]*r
    inverse_list = [1]*r

    for i in xrange(1, r):
        partial_products[i] = partial_products[i - 1]*moduli[i - 1]
        inverse_list[i] = inverse_mod(partial_products[i], moduli[i])

    product = partial_products[r - 1]*moduli[r - 1]
    data = (r, partial_products, inverse_list, product)
    return data


def chinese_preconditioned(L, preconditioning_data=None):
    """
    Returns a solution to the congruences given in L.

    Given a list of congruenes as (residue, modulus) pairs, where all moduli are
    pairwise coprime, this returns the unique solution n in [0, M) to the
    congruences, where M is the product of the moduli.

    Input:
        * L: list
            L is a list of congruences of the form (residue, modulus).

        * preconditioning_data: tuple (default=None)
            This is a tuple containing the preconditioning data for this set of
            moduli. This should be the output of the function
            crt_preconditioning_data.
    
    Output:
        * x: int
            x is the least positive integer satisfying all of the congruences.

    Examples:
        >>> chinese([(0, 4), (1, 25)])
        76

    Details:
        This algorithm is due to Garner. See Algorithm 2.1.7 in "Prime Numbers -
        A Computational Perspective" by Crandall and Pomerance for details.

        If multiple systems of congruences are to be solved, where each system
        has the same set of moduli as the other, then this method is much faster
        than the function "chinese", as long as the preconditioning data is only
        computed once.
        
        A typical use of this function is in combining solutions of some
        polynomial congruence modulo prime powers to find the solutions modulo
        the product of these prime-power moduli. In this case, you will have a
        fixed set of moduli, but multiple residues for each modulus.
    """
    moduli = [m for (a, m) in L]
    residues = [a for (a, m) in L]
    x = residues[0]

    (r, partial_products, inverse_list, product) = preconditioning_data
    for i in xrange(1, r):
        u = (residues[i] - x)*inverse_list[i] % moduli[i]
        x = x + u*partial_products[i]
    return x % product


def jacobi_symbol(a, m):
    """
    Returns the Jacobi symbol (a|m).

    Given positive integers a and m with m odd, this algorithm returns the
    Jacobi symbol (a|m), which for m an odd prime is also the Legendre symbol.

    Input:
        * a: int (a > 0)
        * m: int (m > 0) (odd)

    Output:
        * t: int
            The integer t = (a|m).

    Examples:
        >>> jacobi_symbol(6477, 11213)
        1
        >>> 707*707 % 11213
        6477
        >>> jacobi_symbol(3, 101)
        -1
        >>> jacobi_symbol(7, 77)
        0

    Details:
        For m an odd prime, this is the Legendre symbol. In this case, we have
        (a|m) = 1 if a is a quadratic residue (mod m), -1 if a is a quadratic
        nonresidue (mod m), and 0 if m divides a.

        For m odd, but not necessarily prime, this is the Jacobi symbol. In this
        case, (a|m) = -1 means that a is a quadratic nonresidue (mod m).

        The algorithm we have implemented follows the exposition given in
        Chapter 5 of "Algorithmic Number Theory - Efficient Algorithms" by Bach
        and Shallit. Also, see Chapter 12 of "A Computational Introduction to
        Number Theory and Algebra" by Shoup for details. 
    """
    if m % 2 == 0:
        raise ValueError('m must be odd')
    
    a = a % m
    t = 1
    while a:
        while a % 2 == 0:
            a = a//2
            if m % 8 in (3, 5):
                t = -t
        a, m = m, a
        if a % 4 == 3 and m % 4 == 3:
            t = -t
        a = a % m
    if m == 1:
        return t
    return 0


def valuation(p, n):
    """
    Returns the highest power of p dividing n.

    Input:
        * p: int
        * n: int

    Output:
        * k: int
            The largest integer such that p**k divides n.

    Examples:
        >>> valuation(2, 12)
        2
        >>> valuation(3, 3628800)
        4

    Details:
        This algorithm performs a binary search to find the highest power of p
        that divides n. The method is essentially the same as the one used to
        compute an integer logarithm.
    """
    if n % p != 0 or p > n:
        return 0

    p = abs(p)
    n = abs(n)
    pk = p
    hi = 1
    # Look at b, b^2, b^4, b^8,... to find in which interval a lives
    while n % pk == 0:
        pk = pk**2
        lo = hi
        hi *= 2

    # Now we know that b^lo <= a <= b^hi perform a binary search on this
    # interval to find the exact value n so that b^n <= a < b^(n + 1)
    while hi - lo > 1:
        mid = (lo + hi)//2
        if n % (p**mid) == 0:
            lo = mid
        else:
            hi = mid
    return lo

