# Polynomials class

from rosemary.number_theory.elementary import gcd

class Polynomial(object):
    def __init__(self, L=None, name='x'):
        """
        Initialize a Polynomial object.
        
        This initializes a new Polynomial object with coefficients from L, and
        variable name.

        Input:
            * L - A list or tuple of integers.
            * x - A string.

        Output:
            * f - A Polynomial object.

        Details:
            L can be a list or tuple of numbers, or nothing at all. If L is
            None, then f is initialized to be a polynomial of degree 1 with no
            constant term. If L is [], then f is initialized to be the zero
            polynomial. If len(L) > 0, then f is Initialized to be the
            polynomial with coefficients L[0] + L[1]*x + ... + L[n]*x^n.

            The string "name" only affects how the polynomial is printed. The
            default is "x"
        """
        if isinstance(L, (list, tuple)):
            coeffs = list(L)
            ll = len(coeffs)
            if ll == 0:
                coeffs [0]
                ll = 1
            else:
                while ll > 1 and coeffs[ll - 1] == 0:
                    coeffs.pop()
                    ll -= 1
        else:
            if L is None:
                coeffs = [0, 1]
                ll = 2
            else:
                coeffs = [L]
                ll = 1
        if ll == 0:
            if coeffs[0] == 0:
                self.degree = -1
            else:
                self.degree = 0
        else:
            self.degree = ll - 1
        self._coefficients = coeffs
        self.lead_coeff = coeffs[-1]
        self.variable = name

    def __repr__(self):
        """
        This returns the representation of self for printing purposes.
        """
        if len(self._coefficients) == 1:
            return str(self._coefficients[0])
        terms = [(i, self._coefficients[i]) for i in range(self.degree + 1) if
                self._coefficients[i]]
        terms.sort(reverse = True)

        def print_term(deg, coeff):
            if deg == 0:
                return str(coeff)
            else:
                t_term = []
                if coeff != 1:
                    t_term.append('%s*' % coeff)
                t_term.append('%s' % self.variable)
                if deg > 1:
                    t_term.append('^%s' % deg)
            s = ''.join(t_term)
            return s

        all_terms = []
        if len(terms) == 0:
            return '0'
        elif terms[0][1] < 0:
            all_terms.append('-')
        for (deg, coeff) in terms:
            if deg < self.degree:
                if coeff < 0:
                    all_terms.append(' - ')
                else:
                    all_terms.append(' + ')
            all_terms.append(print_term(deg, abs(coeff)))
        s = ''.join(all_terms)
        return s

    def __add__(self, other):
        """
        Adds the two polynomials.

        This returns the polynomial self + other.
        """
        if not isinstance(other, Polynomial):
            try:
                other_poly = Polynomial(other, self.variable)
                return self + other_poly
            except:
                raise ValueError('Other cannot be coerced into a polynomial')

        coeffs = list_add(self._coefficients, other._coefficients)
        ss = Polynomial(coeffs, self.variable)
        return ss

    def __radd__(self, other):
        """
        Adds the two polynomials.

        This returns the polynomial other + self.
        """
        return self + other

    def __neg__(self):
        """
        Negates self.

        This returns the polynomial -self.
        """
        coeffs = [-e for e in self._coefficients]
        pp = Polynomial(coeffs, self.variable)
        return pp

    def __sub__(self, other):
        """
        Subtracts two polynomials.

        This returns the polynomial self - other.
        """
        if not isinstance(other, Polynomial):
            try:
                other_poly = Polynomial(other, self.variable)
                return self + other_poly
            except:
                raise ValueError('Other cannot be coerced into a polynomial')

        coeffs = list_sub(self._coefficients, other._coefficients)
        ss = Polynomial(coeffs, self.variable)
        return ss

    def __rsub__(self, other):
        """
        Subtracts two polynomials.

        This returns the polynomial other - self.
        """
        pp = -self + other
        return pp

    def __mul__(self, other):
        """
        Multiplies two polynomials.

        This returns the polynomial self * other. The algorithm used is either
        the gradeschool method, or the karatsuba algorithm, depending on the
        degrees of self and other.
        """
        if not isinstance(other, Polynomial):
            try:
                coeffs = [other * c for c in self._coefficients]
                return Polynomial(coeffs, self.variable)
            except:
                raise ValueError('Other cannot be coerced into a polynomial')

        coeffs = mul_karatsuba(self._coefficients, other._coefficients)
        pp = Polynomial(coeffs, self.variable)
        return pp

    def __rmul__(self, other):
        """
        Multiplies two polynomials.

        This returns the polynomial other * self.
        """
        return self * other

    def __pow__(self, k):
        """
        Exponentiates self.

        This returns the polynomial self**k. The algorithm used is the J.C.P
        Miller Power Recurrence.
        """
        coeffs = power_jcp(self._coefficients, k)
        pp = Polynomial(coeffs, self.variable)
        return pp

    def __getitem__(self, index):
        """
        Returns the ith coefficient of self.
        """
        cc = self._coefficients[index]
        return cc

    def __eq__(self, other):
        """
        Test equality between self and other.

        This compares the coefficients of self and other. The variable name is
        ignored.
        """
        if not isinstance(other, Polynomial):
            try:
                return Polynomial(other, self.variable) == self
            except:
                raise ValueError('Other cannot be coerced into a polynomial')
        bb = (self._coefficients == other._coefficients)
        return bb

    def __ne__(self, other):
        """
        Test inequality between self and other.

        This compares the coefficients of self and other. The variable name is
        ignored.
        """
        if not isinstance(other, Polynomial):
            try:
                return Polynomial(other, self.variable) != self
            except:
                raise ValueError('Other cannot be coerced into a polynomial')
        bb = (self._coefficients != other._coefficients)
        return bb

    def __call__(self, x):
        """
        Returns the value of self at x.

        This uses Horner's Rule to evaluate the Polynomial self at the value x.
        """
        ss = 0
        for c in reversed(self._coefficients):
            ss = ss * x + c
        return ss

    def coefficients(self):
        """
        Returns a copy of the coefficients of self.

        This returns a copy of the coefficients of the polynomial. If direct
        access the the coefficients is desired, use the self._coefficients
        attribute.
        """
        coeffs = list(self._coefficients)
        return coeffs

    def content(self):
        """
        Returns the content of self.

        The content of a polynomial with integral coefficients is defined to be
        the greatest common divisor of the coefficients.
        """
        coeffs = list(self._coefficients)
        return gcd(coeffs)

    def primitive_part(self):
        """
        Returns the primitive part of self.

        The primitive part of a polynomial with integral coefficients is defined
        to be the polynomial obtained by dividing out the content of self. It is
        conventional to define the primitive part so that its leading
        coefficient is positive.
        """
        coeffs = list(self._coefficients)
        cc = gcd(coeffs)
        if coeffs[-1] < 0:
            cc *= -1
        pp_coeffs = [ e / cc for e in coeffs ]
        pp = Polynomial(pp_coeffs, self.variable)
        return pp

    def power_binary(self, k):
        """
        Exponentiates self.

        This return self**k. The algorithm used is a binary powering algorithm.
        This is slower than the J.C.P. Miller recurrence, but is included as a
        reference.
        """
        z = self
        y = Polynomial([1], self.variable)
        if k == 0:
            return y
        while True:
            if k % 2 == 1:
                y = z * y
            k = k // 2
            if k == 0:
                return y
            z = z * z

    def pseudo_quo_rem(self, other):
        """
        Performs pseudo division on self and other

        Given polynomials u(x) (self) and v(x) (other), where deg(self) >=
        deg(other) >= 0, this algorithm finds polynomials q(z) and r(x)
        satisfying l**(m - n + 1) * u(x) = q(x) * v(x) + r(x), where l is the
        leading coefficient of v(x).

        Input:
            * self - A Polynomial of degree m >= 0.
            * other - A Polynomial of degree 0 <= n <= m.

        Output:
            * (quo, rem) - quotient and remainder of the division of self by
              other.

        Details:
            The algorithm used is given in Section 4.6.1 of 
        """
        if not isinstance(other, Polynomial):
            try:
                other_poly = Polynomial(other, self.variable)
                return self.pseudo_quo_rem(other_poly)
            except:
                raise ValueError('Other cannot be coerced into a polynomial')
        u = list(self._coefficients)
        v = list(other._coefficients)
        m = len(u) - 1
        n = len(v) - 1
        if m < n:
            qq = Polynomial([0], self.variable)
            return (qq, self)
        q = [0] * (m - n + 1)
        for k in xrange(m - n, -1, -1):
            q[k] = u[n + k] * v[n]**k
            for j in xrange(n + k - 1, -1, -1):
                if j < k:
                    u[j] = v[n] * u[j]
                else:
                    u[j] = v[n] * u[j] - u[n + k] * v[j - k]
        qq = Polynomial(q, self.variable)
        rr = Polynomial(u[0:n], self.variable)
        return (qq, rr)

    def quo_rem(self, other):
        """
        Returns the quotient and remainder of the division.
        
        Given polynomials u(x) (self) and v(x) (other), where deg(self) >=
        deg(other) >= 0, this algorithm finds polynomials q(z) and r(x)
        satisfying u(x) = q(x) * v(x) + r(x).
        """
        if not isinstance(other, Polynomial):
            try:
                other_poly = Polynomial(other, self.variable)
                qq = self.quo_rem(other_poly)
                return qq
            except:
                raise ValueError('Other cannot be coerced into a polynomial')
        a = list(self._coefficients)
        b = list(other._coefficients)
        k = len(a)
        l = len(b)
        if k < l:
            qq = Polynomial([0], self.variable)
            return (qq, g)
        lc = b[l - 1]
        r = a[0:k]
        q = [0] * (k - l + 1)
        if a[k - 1] % lc != 0:
            raise ValueError("Non-Integral Coefficients")
        for i in xrange(k - l, -1, -1):
            q[i] = r[i + l - 1] // lc
            for j in xrange(l):
                r[i + j] -= q[i] * b[j]
        qq = Polynomial(q, self.variable)
        rr = Polynomial(r, self.variable)
        return (qq, rr)

    def gcd(self, other):
        """
        Returns the greatest common divisor of self and other.

        Given polynomials self and other, this returns their greatest common
        divisor, using only exact integral arithmetic of the coefficients. The
        algorithm used is the Subresultant GCD algorithm. For details, see
        Algorithm 3.3.1 of "A Course in Computational Algebraic Number Theory"
        by Cohen, or Algorithm C in section 4.6.1 of "The Art of Computer
        Programmming: Vol II" by Knuth.
        """
        d = gcd(self.content(), other.content())
        u = self.primitive_part()
        v = other.primitive_part()
        h = g = 1
        dd = u.degree - v.degree

        while True:
            (q, r, e) = u.pseudo_quo_rem(v)
            if r == 0:
                return d * v.primitive_part()
            elif r.degree == 0:
                return Polynomial([d])
            u = v
            v = Polynomial([e // (g*h**dd) for e in r._coefficients])
            g = u.lead_coeff
            if dd > 1:
                h = g**dd // h**(dd - 1)
            else:
                h = h**(1 - dd) // g**(-dd)

def mul_classical(M, N):
    """
    Returns the product of M and N.

    Given lists M and N of numbers, this returns the coefficients of the product
    f(x) * g(x), where f(x) = M[0] + M[1]*x + ... + M[k - 1]*x^{k - 1} and g(x)
    = N[0] + N[1]*x + ... + N[l - 1]*x^{l - 1}.
    
    Input:
        * M - A list or tuple of numbers.
        * N - A list or tuple of numbers.

    Output:
        * L - A list of numbers. These are the coefficients of the product of
        f(x) * g(x), where M and N are the lists of their coefficients,
        respectively.

    Details:
        This function uses the classical "gradeschool" algorithm for
        multiplication. This algorithm requires O(k*l) operations to compute the
        product of two polynomials with degrees k >= 1 and l >= 1.

    Examples:
        >>> M = [1, 1]
        >>> mul_classical(M, M)
        [1, 2, 1]
    """
    m = len(M)
    n = len(N)
    if m == 0 or n == 0:
        return []
    P = [0] * (m + n - 1)
    for i in xrange(m):
        for j in xrange(n):
            P[i + j] += M[i] * N[j]
    return P

def mul_karatsuba(M, N):
    """
    Returns the product of M and N.

    Given lists M and N of numbers, this returns the coefficients of the product
    f(x) * g(x), where f(x) = M[0] + M[1]*x + ... + M[n]*x^n and g(x) = N[0] +
    N[1]*x + ... + N[k]*x^k.
     
    Input:
        * M - A list or tuple of numbers.
        * N - A list or tuple of numbers.

    Output:
        * L - A list of numbers. These are the coefficients of the product of
        f(x) * g(x), where M and N are the lists of their coefficients,
        respectively.

    Details:
        This function uses the Karatsuba algorithm for multiplication. For
        multiplying two degree n polynomials, this method has time complexity
        O(n^c), where c = log_2(3) = 1.58496..., so it is much faster than the
        classical O(n^2) method. For more information, see section 4.3.3 of "The
        Art of Computer Programming, Volume II" by D. Knuth.

    Examples:
        >>> M = [1] * 2000
        >>> N = [1] * 3000
        >>> time A = polynomials.mul_classical(M, N)
        CPU times: user 1.18 s, sys: 0.00 s, total: 1.18 s
        Wall time: 1.17 s
        >>> time B = polynomials.mul_karatsuba(M, N)
        CPU times: user 0.30 s, sys: 0.00 s, total: 0.30 s
        Wall time: 0.30 s
        >>> A == B
        True
    """
    m = len(M)
    n = len(N)
    kthresh = 30

    # Classical Multiplication is faster until a certain point.
    if m < kthresh or n < kthresh:
        return mul_classical(M, N)

    # The algorithm works best when we split into halves
    k = max(m, n) // 2

    # write M = f1*x^k + f0, N = g1*x^k + g0
    f0 = M[0:k]
    f1 = M[k:]
    g0 = N[0:k]
    g1 = N[k:]

    # We want to compute MN = (f1*x^k + f0)(g1*x^k + g0).
    # We let z0 = f0*g0, z2 = f1*g1
    z0 = mul_karatsuba(f0, g0)
    z2 = mul_karatsuba(f1, g1)
    f1_f0 = list_add(f1, f0)
    g1_g0 = list_add(g1, g0)
    xxyy = mul_karatsuba(f1_f0, g1_g0)
    z1 = list_sub(list_sub(xxyy, z2), z0)
    t2 = [0] * (2*k) + z2
    t1 = [0] * (k) + z1
    xy = list_add(list_add(t2, t1), z0)
    return xy

def power_jcp(A, k):
    """
    Returns the kth power

    Given a nonnegative integer k and a list A, this returns the coefficients of
    f(x)^k, where f(x) is the polynomial f(x) = A[0] + A[1]*x + ... + A[n]*x^n.

    Input:
        * A - A list or tuple of numbers
        * k - a nonnegative integer

    Output:
        * P - A list of numbers. These are the coefficients of f(x)^k, where
        f(x) is the polynomial with coefficients given by the list A.

    Details:
        The algorithm used is the J.C.P. Miller Pure-Power Recurrence. For a
        fixed polynomial of degree L, this algorithm has time complexity
        O(L^2*m) = O(m) and space complexity O(L) = O(1). For more details, look
        at the paper titled "The J.C.P. Miller Recurrence for Exponentiating a
        Polynomial, and its q-Analog" by D. Zeilberger.

    Examples:
        >>> f = Polynomial([1, 3, 2], 'x')
        >>> f
        2*x^2 + 3*x + 1
        >>> f**2
        4*x^4 + 12*x^3 + 13*x^2 + 6*x + 1
        >>> L = [1, 3, 2]
        >>> power_jcp(L, 2)
        [1, 6, 13, 12, 4]
    """
    shift = 0
    while A[shift] == 0:
        shift += 1

    B = A[shift:]
    n = len(B)
    new_deg = (n - 1) * k + 1
    P = [0] * new_deg
    a0 = B[0]
    P[0] = a0**k
    for i in xrange(1, new_deg):
        ss = 0
        for j in xrange(1, min(n, i + 1)):
            ss += B[j] * ((k + 1) * j - i) * P[i - j]
        P[i] = ss // (i * a0)
    return [0]*shift*k + P

def list_add(A, B):
    """
    Adds two lists.

    Let A be a list of length m, and let B be a list of length n. This returns
    the list [A[i] + B[i] | 0 <= i <= max(m, n)], where A[i], B[i] = 0 if i >=
    m, n, respectively. In effect, this works the same as polynomial addition,
    where we assume that coefficients not appearing in a polynomial are zero.

    Input:
        * A - A list of numbers.
        * B - A list of numbers.

    Output:
        * C - A list of numbers.

    Examples:
        >>> list_add([0, 1, 2], [1, -1, 1])
        [1, 0, 3]
        >>> list_add([1, 1, 1, 1], [2, 3, 5])
        [3, 4, 6, 1]
    """
    if len(A) >= len(B):
        C = list(A)
        for (i, c) in enumerate(B):
            C[i] += c
    else:
        C = list(B)
        for (i, c) in enumerate(A):
            C[i] += c
    return C

def list_sub(A, B):
    """
    Subtract two lists.

    Let A be a list of length m, and let B be a list of length n. This returns
    the list [A[i] - B[i] | 0 <= i <= max(m, n)], where A[i], B[i] = 0 if i >=
    m, n, respectively. In effect, this works the same as polynomial addition,
    where we assume that coefficients not appearing in a polynomial are zero.

    Input:
        * A - A list of numbers.
        * B - A list of numbers.

    Output:
        * C - A list of numbers.

    Examples:
        >>> list_sub([0, 1, 2], [1, -1, 1])
        [-1, 2, 1]
        >>> list_sub([1, 1, 1, 1], [2, 3, 5])
        [-1, -2, -4, 1]
    """

    if len(A) >= len(B):
        C = list(A)
        for (i, c) in enumerate(B):
            C[i] -= c
    else:
        C = [-e for e in B]
        for (i, c) in enumerate(A):
            C[i] += c
    return C

