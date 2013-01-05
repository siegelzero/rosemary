# Class of Integral Binary Quadratic Forms

from rosemary.number_theory.elementary import gcd, integer_sqrt
import math

class BinaryQuadraticForm:
    """
    Class for Binary Quadratic Forms

    Input:
        * "L" - a list or tuple with three entries [a, b, c] or (a, b, c)

    Output:
        * The binary quadratic form a*x^2 + b*x*y + c*y^2

    Examples:
        >>> B = BinaryQuadraticForm([1, 2, 2])
        >>> B
        x^2 + 2*x*y + 2*y^2
        >>> B.discriminant()
        -4
    """
    def __init__(self, abc):
        if isinstance(abc, (list, tuple)):
            if len(abc) == 3:
                self.a, self.b, self.c = abc
        else:
            raise TypeError("Binary quadratic form must be given by a list or tuple of three coefficients a, b, c.")

    def __repr__(self):
        s = ''
        if self._a != 0:
            if self._a < 0:
                s += '-'
            aa = abs(self._a)
            if aa > 1:
                s += str(aa)
                s += '*'
            s += 'x^2'

        if self._b != 0:
            if self._b > 0:
                s += ' + '
            else:
                s += ' - '
            ab = abs(self._b)
            if ab > 1:
                s += str(ab)
                s += '*'
            s += 'x*y'

        if self._c != 0:
            if self._c > 0:
                s += ' + '
            else:
                s += ' - '
            ac = abs(self._c)
            if ac > 1:
                s += str(ac)
                s += '*'
            s += 'y^2'

        if s != '':
            return s
        else:
            return '0'

    def __call__(self, x, y = None):
        if y is None and isinstance(x, (list, tuple)):
            x, y = x[0], x[1]

        return self._a * x * x + self._b * x * y + self._c * y * y

    def __getitem__(self, i):
        """
        Returns the ith coefficient of self

        Examples:
        >>> B = BinaryQuadraticForm([1, 4, 3])
        >>> B[0]
        1
        """
        return (self.a, self.b, self.c)[i]

    ################################################################################
    # Methods for the class
    ################################################################################

    def discriminant(self):
        """
        Returns the discriminant d = b^2 - 4*a*c of the binary quadratic form
        a*x^2 + b*x*y + c*y^2.
        """
        return self.b**2 - 4 * self.a * self.c

    def is_primitive(self):
        """
        Returns whether or not the binary quadratic form self is primitive: The
        binary quadratic form a*x^2 + b*x*y + c*y^2 is primitive
        if gcd(a, b, c) = 1.
        """
        return number_theory.gcd(list(self)) == 1

    def is_reduced(self):
        """
        Returns whether or not the binary quadratic form is reduced. The BQF
        a*x^2 + b*x*y + c*y^2 is reduced if -a < b <= a < c or 0 <= |b| <= a = c
        """
        return ((-self.a < self.b <= self.a < self.c) or 0 <= abs(self.b) <=
                self.a == self.c)
    
    def representations(self, n):
        a = self.a
        b = self.b
        c = self.c
        d = abs(b*b - 4*a*c)
        L = []
        x = 0
        while x*x <= 4*c*n // d:
            y = 0
            while y*y <= 4*a*n // d:
                tt = a*x*x + c*y*y
                if tt + b*x*y == n:
                    L += [ (x, y), (-x, -y) ]
                elif tt - b*x*y == n:
                    L += [ (x, -y), (-x, y) ]
                y += 1
            x += 1

        return L
 
################################################################################

def class_number(d):
    """
    class_number(d):
    Given a discriminant d < 0, d = 0, 1 (mod 4), this returns the class number
    h(d); i.e. the number of proper equivalence classes of binary quadratic forms
    of discriminant d.
    """
    if d >= 0:
        raise ValueError("discriminant d must be negative")
    if d % 4 not in (0, 1):
        raise ValueError("discriminant d must satisfy d = 0, 1 (mod 4)")

    h = 0
    b = d % 2
    # we want -a < b <= a < c or 0 <= b <= a = c
    while b*b <= abs(d) // 3:
        A = (b*b - d) // 4
        a = max(1, b)
        while a*a <= A:
            if A % a == 0:
                c = A // a
                if number_theory.gcd([a, b, c]) == 1:
                    if b == 0 or a in (b, c):
                        h += 1
                    else:
                        h += 2
            a += 1
        b += 2

    return h

################################################################################

def cornacchia(d, p):
    """
    cornacchia(d, p):
    Given a prime number p and an integer d such that 0 < d < p, this algorithm
    either outputs an integer solution (x, y) to the Diophantine equation
    x^2 + d*y^2 = p, or says that such a solution does not exist.
    """
    k = number_theory.jacobi_symbol(-d, p)
    if k == -1:
        raise ValueError("The equation has no solutions")

    x0 = number_theory.sqrt_mod_p(-d, p)

    if x0 <= p // 2:
        x0 = p - x0

    a = p
    b = x0
    l = number_theory.integer_sqrt(p)

    while b > l:
        r = a % b
        a = b
        b = r

    t = p - b*b

    if t % d != 0:
        raise ValueError("The equation has no solutions")

    c = t // d
    s = number_theory.integer_sqrt(c)
    if s*s != c:
        raise ValueError("The equation has no solutions")

    return (b, s)

################################################################################
