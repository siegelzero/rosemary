import rosemary.number_theory.sieves

from rosemary.number_theory.core import integer_nth_root, integer_sqrt
from rosemary.data_structures import bit_sieve
from rosemary.number_theory.prime_list import _PRIME_LIST

from bisect import bisect
from collections import defaultdict


###########################################################################
# Algorithms for counting primes.
###########################################################################

def legendre(n):
    r"""Returns the number of primes <= n.

    Parameters
    ----------
    n : int (n > 0)

    Returns
    -------
    count : int
        The number of primes <= n.

    Notes
    -----
    This method uses the function `primes` to compute the primes, and
    simple does some bookkeeping to maintain a running count.

    Examples
    --------
    >>> pi_table(10)
    ([2, 3, 5, 7], [0, 0, 1, 2, 2, 3, 3, 4, 4, 4, 4])
    """
    if n <= 1:
        return 0
    elif n <= _PRIME_LIST[-1]:
        return bisect(_PRIME_LIST, n)

    root = integer_sqrt(n)
    primes = rosemary.number_theory.sieves.primes(root)
    a = len(primes)

    def phi(x, a, cache={}):
        if (x, a) in cache:
            return cache[x, a]

        if a == 1:
            return (x + 1)//2

        value = (x + 1)//2
        for i in xrange(2, a + 1):
            p = primes[i - 1]
            if p > x:
                break
            value -= phi(x//p, i - 1)

        cache[x, a] = value
        return value

    return phi(n, a) + a - 1


def meissel(n):
    if n <= 1:
        return 0
    elif n <= _PRIME_LIST[-1]:
        return bisect(_PRIME_LIST, n)

    root = integer_nth_root(3, n**2)
    primes, pi_cache = pi_table(root)

    a = pi_cache[integer_nth_root(3, n)]
    b = pi_cache[integer_nth_root(2, n)]

    def phi(x, a, cache={}):
        if (x, a) in cache:
            return cache[x, a]

        b = pi_cache[int((x + 0.5)**(0.5))]
        c = pi_cache[int((x + 0.5)**(1.0/3.0))]

        if a <= 7:
            if a == 1:
                value = x - x//2

            elif a == 2:
                value = x - x//2 - x//3 + x//6

            elif a == 3:
                value = x - x//2 - x//3 - x//5 + x//6 + x//10 + x//15 - x//30

            elif a == 4:
                value = (x - x//2 - x//3 - x//5 - x//7 + x//6 + x//10 + x//14 + x//15 + x//21 + x//35 - x//30 - x//42 -
                         x//70 - x//105 + x//210)
            elif a == 5:
                value = (x - x//2 - x//3 - x//5 - x//7 - x//11 + x//6 + x//10 + x//14 + x//22 + x//15 + x//21 + x//33 +
                         x//35 + x//55 + x//77 - x//30 - x//42 - x//66 - x//70 - x//110 - x//154 - x//105 - x//165 -
                         x//231 - x//385 + x//210 + x//330 + x//462 + x//770 + x//1155 - x//2310)
            elif a == 6:
                value = (x - x//2 - x//3 - x//5 - x//7 - x//11 - x//13 + x//6 + x//10 + x//14 + x//22 + x//26 + x//15 +
                         x//21 + x//33 + x//39 + x//35 + x//55 + x//65 + x//77 + x//91 + x//143 - x//30 - x//42 - x//66
                         - x//78 - x//70 - x//110 - x//130 - x//154 - x//182 - x//286 - x//105 - x//165 - x//195 -
                         x//231 - x//273 - x//429 - x//385 - x//455 - x//715 - x//1001 + x//210 + x//330 + x//390 +
                         x//462 + x//546 + x//858 + x//770 + x//910 + x//1430 + x//2002 + x//1155 + x//1365 + x//2145 +
                         x//3003 + x//5005 - x//2310 - x//2730 - x//4290 - x//6006 - x//10010 - x//15015 + x//30030)
            else:
                value = (x - x//2 - x//3 - x//5 - x//7 - x//11 - x//13 - x//17 + x//6 + x//10 + x//14 + x//22 + x//26 +
                         x//34 + x//15 + x//21 + x//33 + x//39 + x//51 + x//35 + x//55 + x//65 + x//85 + x//77 + x//91 +
                         x//119 + x//143 + x//187 + x//221 - x//30 - x//42 - x//66 - x//78 - x//102 - x//70 - x//110 -
                         x//130 - x//170 - x//154 - x//182 - x//238 - x//286 - x//374 - x//442 - x//105 - x//165 -
                         x//195 - x//255 - x//231 - x//273 - x//357 - x//429 - x//561 - x//663 - x//385 - x//455 -
                         x//595 - x//715 - x//935 - x//1105 - x//1001 - x//1309 - x//1547 - x//2431 + x//210 + x//330 +
                         x//390 + x//510 + x//462 + x//546 + x//714 + x//858 + x//1122 + x//1326 + x//770 + x//910 +
                         x//1190 + x//1430 + x//1870 + x//2210 + x//2002 + x//2618 + x//3094 + x//4862 + x//1155 +
                         x//1365 + x//1785 + x//2145 + x//2805 + x//3315 + x//3003 + x//3927 + x//4641 + x//7293 +
                         x//5005 + x//6545 + x//7735 + x//12155 + x//17017 - x//2310 - x//2730 - x//3570 - x//4290 -
                         x//5610 - x//6630 - x//6006 - x//7854 - x//9282 - x//14586 - x//10010 - x//13090 - x//15470 -
                         x//24310 - x//34034 - x//15015 - x//19635 - x//23205 - x//36465 - x//51051 - x//85085 +
                         x//30030 + x//39270 + x//46410 + x//72930 + x//102102 + x//170170 + x//255255 - x//510510)

        elif b > a >= c and x <= root:
            value = pi_cache[x] - a + 1 + phi2(x, a, b)

        elif a >= b and x <= root:
            value = pi_cache[x] - a + 1

        else:
            value = (x + 1)//2
            for i in xrange(1, a):
                if primes[i] > x:
                    break
                value -= phi(x//primes[i], i)

        cache[x, a] = value
        return value

    def phi2(n, a, b):
        total = (a - b)*(b + a - 1)//2
        for i in xrange(a, b):
            total += pi_cache[n//primes[i]]
        return total

    value = phi(n, a) - 1 + a - phi2(n, a, b)
    return value


def lehmer(n):
    if n <= 1:
        return 0
    elif n <= _PRIME_LIST[-1]:
        return bisect(_PRIME_LIST, n)

    root = integer_nth_root(4, n**3)
    primes = rosemary.number_theory.sieves.primes(root)

    a = bisect(primes, integer_nth_root(4, n))
    b = bisect(primes, integer_sqrt(n))
    c = bisect(primes, integer_nth_root(3, n))

    def phi(x, a, cache={}):
        if (x, a) in cache:
            return cache[x, a]

        if a <= 7:
            if a == 1:
                value = x - x//2

            elif a == 2:
                value = x - x//2 - x//3 + x//6

            elif a == 3:
                value = x - x//2 - x//3 - x//5 + x//6 + x//10 + x//15 - x//30

            elif a == 4:
                value = (x - x//2 - x//3 - x//5 - x//7 + x//6 + x//10 + x//14 + x//15 + x//21 + x//35 - x//30 - x//42 -
                         x//70 - x//105 + x//210)
            elif a == 5:
                value = (x - x//2 - x//3 - x//5 - x//7 - x//11 + x//6 + x//10 + x//14 + x//22 + x//15 + x//21 + x//33 +
                         x//35 + x//55 + x//77 - x//30 - x//42 - x//66 - x//70 - x//110 - x//154 - x//105 - x//165 -
                         x//231 - x//385 + x//210 + x//330 + x//462 + x//770 + x//1155 - x//2310)
            elif a == 6:
                value = (x - x//2 - x//3 - x//5 - x//7 - x//11 - x//13 + x//6 + x//10 + x//14 + x//22 + x//26 + x//15 +
                         x//21 + x//33 + x//39 + x//35 + x//55 + x//65 + x//77 + x//91 + x//143 - x//30 - x//42 - x//66
                         - x//78 - x//70 - x//110 - x//130 - x//154 - x//182 - x//286 - x//105 - x//165 - x//195 -
                         x//231 - x//273 - x//429 - x//385 - x//455 - x//715 - x//1001 + x//210 + x//330 + x//390 +
                         x//462 + x//546 + x//858 + x//770 + x//910 + x//1430 + x//2002 + x//1155 + x//1365 + x//2145 +
                         x//3003 + x//5005 - x//2310 - x//2730 - x//4290 - x//6006 - x//10010 - x//15015 + x//30030)
            else:
                value = (x - x//2 - x//3 - x//5 - x//7 - x//11 - x//13 - x//17 + x//6 + x//10 + x//14 + x//22 + x//26 +
                         x//34 + x//15 + x//21 + x//33 + x//39 + x//51 + x//35 + x//55 + x//65 + x//85 + x//77 + x//91 +
                         x//119 + x//143 + x//187 + x//221 - x//30 - x//42 - x//66 - x//78 - x//102 - x//70 - x//110 -
                         x//130 - x//170 - x//154 - x//182 - x//238 - x//286 - x//374 - x//442 - x//105 - x//165 -
                         x//195 - x//255 - x//231 - x//273 - x//357 - x//429 - x//561 - x//663 - x//385 - x//455 -
                         x//595 - x//715 - x//935 - x//1105 - x//1001 - x//1309 - x//1547 - x//2431 + x//210 + x//330 +
                         x//390 + x//510 + x//462 + x//546 + x//714 + x//858 + x//1122 + x//1326 + x//770 + x//910 +
                         x//1190 + x//1430 + x//1870 + x//2210 + x//2002 + x//2618 + x//3094 + x//4862 + x//1155 +
                         x//1365 + x//1785 + x//2145 + x//2805 + x//3315 + x//3003 + x//3927 + x//4641 + x//7293 +
                         x//5005 + x//6545 + x//7735 + x//12155 + x//17017 - x//2310 - x//2730 - x//3570 - x//4290 -
                         x//5610 - x//6630 - x//6006 - x//7854 - x//9282 - x//14586 - x//10010 - x//13090 - x//15470 -
                         x//24310 - x//34034 - x//15015 - x//19635 - x//23205 - x//36465 - x//51051 - x//85085 +
                         x//30030 + x//39270 + x//46410 + x//72930 + x//102102 + x//170170 + x//255255 - x//510510)

        elif a >= bisect(primes, x**(0.5)) and x < primes[-1]:
            pi = bisect(primes, x)
            value = pi - a + 1

        else:
            value = (x + 1)//2
            for i in xrange(1, a):
                if primes[i] > x:
                    break
                value -= phi(x//primes[i], i)

        cache[x, a] = value
        return value

    def phi2(n, a, b):
        s1 = 0
        for i in xrange(a, b):
            s1 += phi(n//primes[i], b) + b - 1
        return s1

    def phi3(n, a, b, c):
        s2 = 0
        for i in xrange(a, c):
            bi = phi(integer_sqrt(n//primes[i]), b) + b - 1
            for j in xrange(i, bi):
                idx = phi(n//(primes[i]*primes[j]), b) + b - 1
                s2 += idx - j
        return s2

    value = (b + a - 2)*(b - a + 1)//2 + phi(n, a) - phi2(n, a, b) - phi3(n, a, b, c)
    return value


def lmo(x):
    if x <= 1:
        return 0
    elif x <= _PRIME_LIST[-1]:
        return bisect(_PRIME_LIST, x)

    root = integer_nth_root(3, x**2)
    primes = rosemary.number_theory.sieves.primes(root)

    t = integer_nth_root(3, x)
    c = bisect(primes, t)
    b = bisect(primes, integer_nth_root(2, x))

    value = (b + c - 2)*(b - c + 1)//2

    for i in xrange(c, b):
        idx = bisect(primes, x//primes[i])
        value -= idx

    special = defaultdict(list)

    stack = [(1, c, 1)]
    push = stack.append
    pop = stack.pop

    # print "recursing"

    while stack:
        (n, a, sign) = pop()

        if a == 1 and n <= t:
            if sign > 0:
                value += (x//n + 1)//2
            else:
                value -= (x//n + 1)//2
        elif n > t:
            special[a].append((x//n, sign))
        else:
            push((n, a - 1, sign))
            push((n*primes[a - 1], a - 1, -sign))

    block = [1]*(root + 1)
    block[0] = 0

    for a in sorted(special):
        p = primes[a - 1]

        block[p::p] = [0]*(root//p)
        last_v = 0
        last_s = 0

        for v, sign in sorted(special[a]):
            last_s += sum(block[last_v:v + 1])
            last_v = v + 1

            if sign > 0:
                value += last_s
            else:
                value -= last_s

    return value


def lmo_bit(x):
    root = int(x**(2.0/3.0))

    primes = rosemary.number_theory.sieves.primes(root)
    t = x**(0.33333333333333)

    c = bisect(primes, t)
    b = bisect(primes, x**(0.5))

    value = (b + c - 2)*(b - c + 1)//2

    for i in xrange(c, b):
        idx = bisect(primes, x//primes[i])
        value -= idx

    special = defaultdict(list)

    stack = [(1, c, 1)]
    push = stack.append
    pop = stack.pop

    while stack:
        (n, a, sign) = pop()

        if a == 1 and n <= t:
            if sign > 0:
                value += (x//n + 1)//2
            else:
                value -= (x//n + 1)//2
        elif n > t:
            special[a].append((x//n, sign))
        else:
            push((n, a - 1, sign))
            push((n*primes[a - 1], a - 1, -sign))

    block = bit_sieve.BITSieveArray(root)
    mark = block.mark_multiples
    total = block.partial_sum

    for a in sorted(special):
        p = primes[a - 1]

        mark(p)

        for v, sign in sorted(special[a]):
            if sign > 0:
                value += total(v)
            else:
                value -= total(v)

    return value


###########################################################################
# Algorithms for summing primes.
###########################################################################


def prime_sum2(n):
    root = integer_sqrt(n)
    primes = rosemary.number_theory.sieves.primes(root)
    a = len(primes)

    sum_table = [0]*(root + 1)
    accumulation = primes[:]
    for i in xrange(1, len(accumulation)):
        accumulation[i] += accumulation[i - 1]

    for (p, total) in zip(primes, accumulation):
        sum_table[p] = total

    for i in xrange(2, root + 1):
        if not sum_table[i]:
            sum_table[i] = sum_table[i - 1]

    def phi(x, a, cache={}):
        if (x, a) in cache:
            return cache[x, a]

        if x == 0:
            return 0
        if a == 0:
            return x*(x + 1)//2

        value = x*(x + 1)//2
        for i in xrange(1, a + 1):
            p = primes[i - 1]
            if p > x:
                break
            value -= p*phi(x//p, i - 1)

        cache[x, a] = value
        return value

    return phi(n, a) + sum(primes) - 1


def prime_sum(n):
    r = integer_sqrt(n)
    V = [n//i for i in range(1, r + 1)]
    V += range(V[-1] - 1, 0, -1)
    S = {i: i*(i + 1)//2 - 1 for i in V}

    for p in range(2, r + 1):
        if S[p] > S[p-1]:
            sp = S[p-1]
            p2 = p*p
            for v in V:
                if v < p2:
                    break
                S[v] -= p*(S[v//p] - sp)

    return S[n]


###########################################################################
# Algorithms for generating tables.
###########################################################################


def pi_table(n):
    r"""Returns the primes <= n, as well as a table of pi(x) for x <= n.

    Parameters
    ----------
    n : int

    Returns
    -------
    (primes, table) : tuple
        The list `primes` contains the primes <= n. The list `table` is
        indexed 0 to n, with `table[k]` giving the number of primes <= k.

    Notes
    -----
    This method uses the function `primes` to compute the primes, and
    simple does some bookkeeping to maintain a running count.

    Examples
    --------
    >>> pi_table(10)
    ([2, 3, 5, 7], [0, 0, 1, 2, 2, 3, 3, 4, 4, 4, 4])
    """
    primes = rosemary.number_theory.sieves.primes(n)
    table = [0]*(n + 1)
    last = 0

    for (i, p) in enumerate(primes):
        table[last:p] = [i]*(p - last)
        last = p

    table[last:n + 1] = [i + 1]*(n + 1 - last)

    return primes, table


def prime_count(n):
    r = integer_sqrt(n)
    V = [n//i for i in range(1, r + 1)]
    V += range(V[-1] - 1, 0, -1)
    S = {i: (i + 1)//2 for i in V}

    for p in range(2, r + 1):
        if S[p] > S[p-1]:
            sp = S[p-1]
            p2 = p*p
            for v in V:
                if v < p2:
                    break
                S[v] -= (S[v//p] - sp)
    return S
