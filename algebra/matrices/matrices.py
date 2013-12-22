class Matrix(object):
    def __init__(self, m, n=None):
        """
        Initialize a Matrix object.

        Input signatures:
            * (m, n): (int, int)
                In this case, an m x n zero matrix is returned.

            * m: int
                In this case, an m x m zero matrix is returned.

            * m: list
                In this case, a matrix is created containing the entries of m.

        Output:
            * M: (Matrix)
                A Matrix object.

        Examples:
            >>> Matrix(2, 3)
            [0 0 0]
            [0 0 0]
            >>> Matrix(4)
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]
            >>> Matrix([[2, 3], [5, 7]])
            [2 3]
            [5 7]
        """
        if n is None:
            if isinstance(m, (int, long)):
                entries = [[0]*m for _ in xrange(m)]
            else:
                entries = list(m)
        elif isinstance(n, (int, long)):
            entries = [[0]*n for _ in xrange(m)]

        self.__entries = entries
        self.__nrows = len(entries)
        self.__ncols = len(entries[0])

    def __repr__(self):
        """
        This returns the representation of self.
        """
        max_len = 0
        nrows = self.__nrows
        ncols = self.__ncols
        entries = self.__entries

        # Find the length of the longest element in self for aligment purposes.
        for i in xrange(nrows):
            for j in xrange(ncols):
                item = str(entries[i][j])
                if len(item) > max_len:
                    max_len = len(item)

        items = []
        for i in xrange(nrows):
            items.append('[')
            for j in xrange(ncols):
                element = entries[i][j]
                str_form = str(element).rjust(max_len)
                items.append(str_form)
                if j < ncols - 1:
                    items.append(' ')
            items.append(']')
            if i < nrows - 1:
                items.append('\n')
        rep = ''.join(items)
        return rep

    def __getitem__(self, key):
        if isinstance(key, (int, long)):
            return self.__entries[key]
        elif isinstance(key, (tuple, list)):
            return self.__entries[key[0]][key[1]]

    def __setitem__(self, key, value):
        if isinstance(key, (int, long)):
            self.__entries[key] = value
        elif isinstance(key, (tuple, list)):
            self.__entries[key[0]][key[1]] = value

    def __add__(self, other):
        """
        Adds the two matrices.

        This returns the matrix self + other.
        """

        if not isinstance(other, Matrix):
            try:
                parent_class = self.__class__
                other_poly = parent_class(other, self.variable)
                return self + other_poly
            except:
                raise ValueError('Other cannot be coerced into a polynomial')

        self_ncols = self.__ncols
        self_nrows = self.__nrows
        other_ncols = self.__ncols
        other_nrows = self.__nrows

        if self_ncols != other_ncols or self_nrows != other_nrows:
            raise ValueError("Invalid dimensions for matrix addition")

        entries = self.__entries
        other_entries = other.__entries
        result = [[0]*self_ncols for _ in xrange(self_nrows)]

        for i in xrange(self_nrows):
            for j in xrange(self_ncols):
                result[i][j] = entries[i][j] + other_entries[i][j]

        parent_class = self.__class__
        return parent_class(result)

    def __radd__(self, other):
        """
        Adds the two matrices.

        This returns the matrix other + self.
        """
        result = self + other
        return result

    def __neg__(self):
        """
        Negates self.

        This returns the matrix -self.
        """
        ncols = self.__ncols
        nrows = self.__nrows
        result = [[0]*ncols for _ in xrange(nrows)]
        entries = self.__entries

        for i in xrange(nrows):
            for j in xrange(ncols):
                result[i][j] = -entries[i][j]

        parent_class = self.__class__
        return parent_class(result)

    def __sub__(self, other):
        """
        Subtracts two matrices.

        This returns the matrix self - other.
        """
        parent_class = self.__class__

        if not isinstance(other, Matrix):
            try:
                other_poly = parent_class(other, self.variable)
                return self + other_poly
            except:
                raise ValueError("Other cannot be coerced into a Matrix")

        ncols = self.__ncols
        nrows = self.__nrows
        entries = list(self.__entries)
        other_entries = list(other.__entries)
        result = [[0]*ncols for _ in xrange(nrows)]

        for i in xrange(nrows):
            for j in xrange(ncols):
                result[i][j] = entries[i][j] - other_entries[i][j]

        return parent_class(result)

    def __rsub__(self, other):
        """
        Subtracts two matrices.

        This returns the matrix other - self.
        """
        return -self + other

    def __mul__(self, other):
        """
        Multiplies two matrices.
        """
        parent_class = self.__class__
        ncols = self.__ncols
        nrows = self.__nrows
        entries = self.__entries

        if not isinstance(other, Matrix):
            try:
                result = [[0]*ncols for _ in xrange(nrows)]
                for i in xrange(nrows):
                    for j in xrange(ncols):
                        result[i][j] = other*entries[i][j]
                return parent_class(result)
            except:
                raise ValueError("Other cannot be coerced into a Matrix")

        other_nrows = other.__nrows
        other_ncols = other.__ncols

        if ncols != other_nrows:
            raise ValueError("Incompatible Matrix dimensions")
        
        other_entries = other.__entries
        result = [[0]*nrows for _ in xrange(other_ncols)]

        for i in xrange(nrows):
            for j in xrange(other_ncols):
                total = 0
                for k in xrange(ncols):
                    total += entries[i][k]*other_entries[k][j]
                result[i][j] = total

        return parent_class(result)

    def __rmul__(self, other):
        """
        Multiplies two polynomials.

        This returns the polynomial other*self.
        """
        parent_class = self.__class__
        ncols = self.__ncols
        nrows = self.__nrows
        result = [[0]*ncols for _ in xrange(nrows)]
        entries = self.__entries

        if not isinstance(other, Matrix):
            try:
                for i in xrange(nrows):
                    for j in xrange(ncols):
                        result[i][j] = other*entries[i][j]
                return parent_class(result)
            except:
                raise ValueError("Other cannot be coerced into a Matrix")

        other_nrows = other.__nrows
        other_ncols = other.__ncols

        if other_ncols != nrows:
            raise ValueError("Incompatible Matrix dimensions")
        
        other_entries = other.__entries
        result = [[0]*other_ncols for _ in xrange(nrows)]

        for i in xrange(other_nrows):
            for j in xrange(ncols):
                total = 0
                for k in xrange(other_ncols):
                    total += entries[i][k]*other_entries[k][j]
                result[i][j] = total

        return parent_class(result)

    def delete_col(self, idx):
        result = []

        for row in self.__entries:
            new_row = []
            for (i, entry) in enumerate(row):
                if i == idx:
                    continue
                new_row.append(entry)
            result.append(new_row)

        parent_class = self.__class__
        return parent_class(result)

    def delete_row(self, idx):
        result = []

        for (i, row) in enumerate(self.__entries):
            if i == idx:
                continue
            result.append(row)

        parent_class = self.__class__
        return parent_class(result)

    def determinant(self):
        """
        Returns the determinant of self.

        Given an n x n matrix M with coefficients in an integral domain R, this
        algorithm computes the determinant of M.

        Input:
            * self (Matrix)

        Output:
            * det (same type as entries of self)

        Examples:
            >>>

        Details:
            The algorithm used here is the Gauss-Bareiss algorithm outlined as
            Algorithm 2.2.6 in 'A Course in Computational Algebraic Number
            Theory' by Cohen. One advantage of this algorithm is that only exact
            divisions are performed, so all intermediate results stay the same
            type as the entries of the input matrix.
        """
        den = 1
        sign = 1
        ncols = self.__ncols
        nrows = self.__nrows

        if nrows != ncols:
            raise ValueError("Determinant only defined for square matrices.")

        m = self.__entries

        for k in xrange(ncols):
            pivot = m[k][k]
            if pivot == 0:
                # Find the first nonzero entry m[i][k] in column k.
                for i in xrange(k + 1, ncols):
                    if m[i][k] != 0:
                        break

                # If no nonzero entry was found in column k, then the
                # determinant is zero.
                if i == ncols - 1:
                    return 0
                # Othwerwise, we perform exchanges.
                for j in xrange(k, ncols):
                    m[i][j], m[k][j] = m[k][j], m[i][j]
                sign *= -1
                pivot = m[k][k]
            for i in xrange(k + 1, ncols):
                for j in xrange(k + 1, ncols):
                    # This is exact division.
                    m[i][j] = (pivot*m[i][j] - m[i][k]*m[k][j])//den
            den = pivot
        return sign*m[ncols - 1][ncols - 1]


class MatrixZZ(Matrix):
    pass

class MatrixQQ(Matrix):
    pass

