# Polynomials

from rosemary.number_theory.core import gcd, gcd_list

from rosemary.algebra.polynomials.algorithms import (
    list_add,
    list_subtract,
    mul_karatsuba,
    power_jcp,
)


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

            The string `name` only affects how the polynomial is printed. The
            default is 'x'
        """
        if isinstance(L, (list, tuple)):
            coeffs = list(L)
            num_coeffs = len(coeffs)
            if num_coeffs == 0:
                coeffs = [0]
                num_coeffs = 1
            else:
                while num_coeffs > 1 and coeffs[num_coeffs - 1] == 0:
                    coeffs.pop()
                    num_coeffs -= 1
        else:
            if L is None:
                coeffs = [0, 1]
                num_coeffs = 2
            else:
                coeffs = [L]
                num_coeffs = 1
        if num_coeffs == 0:
            if coeffs[0] == 0:
                self.degree = -1
            else:
                self.degree = 0
        else:
            self.degree = num_coeffs - 1
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
        terms.sort(reverse=True)

        def print_term(deg, coeff):
            if deg == 0:
                return str(coeff)
            else:
                term_pieces = []
                if coeff != 1:
                    term_pieces.append('%s*' % coeff)
                term_pieces.append('%s' % self.variable)
                if deg > 1:
                    term_pieces.append('^%s' % deg)
            term = ''.join(term_pieces)
            return term

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
        result = ''.join(all_terms)
        return result

    def __add__(self, other):
        """
        Adds the two polynomials.

        This returns the polynomial self + other.
        """
        if not isinstance(other, Polynomial):
            try:
                parent_class = self.__class__
                other_poly = parent_class(other, self.variable)
                return self + other_poly
            except:
                raise ValueError('Other cannot be coerced into a polynomial')

        parent_class = self.__class__
        coeffs = list_add(self._coefficients, other._coefficients)
        result = parent_class(coeffs, self.variable)
        return result

    def __radd__(self, other):
        """
        Adds the two polynomials.

        This returns the polynomial other + self.
        """
        result = self + other
        return result

    def __neg__(self):
        """
        Negates self.

        This returns the polynomial -self.
        """
        coeffs = [-e for e in self._coefficients]
        parent_class = self.__class__
        result = parent_class(coeffs, self.variable)
        return result

    def __sub__(self, other):
        """
        Subtracts two polynomials.

        This returns the polynomial self - other.
        """
        if not isinstance(other, Polynomial):
            try:
                parent_class = self.__class__
                other_poly = parent_class(other, self.variable)
                return self + other_poly
            except:
                raise ValueError('Other cannot be coerced into a polynomial')

        parent_class = self.__class__
        coeffs = list_subtract(self._coefficients, other._coefficients)
        result = parent_class(coeffs, self.variable)
        return result

    def __rsub__(self, other):
        """
        Subtracts two polynomials.

        This returns the polynomial other - self.
        """
        result = -self + other
        return result

    def __mul__(self, other):
        """
        Multiplies two polynomials.

        This returns the polynomial self*other. The algorithm used is either the
        gradeschool method or the karatsuba algorithm, depending on the degrees
        of self and other.
        """
        if not isinstance(other, Polynomial):
            try:
                coeffs = [other*c for c in self._coefficients]
                parent_class = self.__class__
                return parent_class(coeffs, self.variable)
            except:
                raise ValueError('Other cannot be coerced into a polynomial')

        parent_class = self.__class__
        coeffs = mul_karatsuba(self._coefficients, other._coefficients)
        result = parent_class(coeffs, self.variable)
        return result

    def __rmul__(self, other):
        """
        Multiplies two polynomials.

        This returns the polynomial other * self.
        """
        result = self*other
        return result

    def __pow__(self, k):
        """
        Exponentiates self.

        This returns the polynomial self**k. The algorithm used is the J.C.P
        Miller Power Recurrence.
        """
        parent_class = self.__class__
        power_coeffs = power_jcp(self._coefficients, k)
        result = parent_class(power_coeffs, self.variable)
        return result

    def __getitem__(self, index):
        """
        Returns the ith coefficient of self.
        """
        coeff = self._coefficients[index]
        return coeff

    def __eq__(self, other):
        """
        Test equality between self and other.

        This compares the coefficients of self and other. The variable name is
        ignored.
        """
        if not isinstance(other, Polynomial):
            try:
                parent_class = self.__class__
                is_equal = (parent_class(other, self.variable) == self)
                return is_equal
            except:
                raise ValueError('Other cannot be coerced into a polynomial')
        is_equal = (self._coefficients == other._coefficients)
        return is_equal

    def __ne__(self, other):
        """
        Test inequality between self and other.

        This compares the coefficients of self and other. The variable name is
        ignored.
        """
        if not isinstance(other, Polynomial):
            try:
                parent_class = self.__class__
                is_equal = (parent_class(other, self.variable) == self)
                return not is_equal
            except:
                raise ValueError('Other cannot be coerced into a polynomial')
        is_equal = (self._coefficients == other._coefficients)
        return not is_equal

    def __call__(self, x):
        """
        Returns the value of self at x.

        This uses Horner's Rule to evaluate the Polynomial self at the value x.
        """
        value = 0
        for coeff in reversed(self._coefficients):
            value = value*x + coeff
        return value

    def coefficients(self):
        """
        Returns a copy of the coefficients of self.

        This returns a copy of the coefficients of the polynomial. If direct
        access the the coefficients is desired, use the self._coefficients
        attribute.
        """
        coeffs = list(self._coefficients)
        return coeffs

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
            g = Polynomial(a, self.variable)
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

        if dd < 0:
            return other.gcd(self)

        while True:
            (q, r) = u.pseudo_quo_rem(v)
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


class PolyZZ(Polynomial):
    def content(self):
        """
        Returns the content of self.

        The content of a polynomial with integral coefficients is defined to be
        the greatest common divisor of the coefficients.
        """
        result = gcd_list(self._coefficients)
        return result

    def primitive_part(self):
        """
        Returns the primitive part of self.

        The primitive part of a polynomial with integral coefficients is defined
        to be the polynomial obtained by dividing out the content of self. It is
        conventional to define the primitive part so that its leading
        coefficient is positive.
        """
        coeffs = self._coefficients
        content = gcd_list(coeffs)
        if coeffs[-1] < 0:
            content *= -1
        parent_class = self.__class__
        primitive_part_coeffs = [coeff//content for coeff in coeffs]
        result = parent_class(primitive_part_coeffs, self.variable)
        return result


class PolyQQ(Polynomial):
    pass
