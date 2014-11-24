import rosemary.number_theory.sieves
from rosemary.number_theory.tables import lookup

from collections import defaultdict
from bisect import bisect

import sys


def legendre(n):
    """
    Returns the number of primes <= n.
    """
    sqrt = int(n**(0.5))
    primes = rosemary.number_theory.sieves.primes(sqrt)
    num = len(primes)

    value = num - 1
    stack = [(1, 0, 1)]
    push = stack.append
    pop = stack.pop

    while stack:
        (prod, k, sign) = pop()
        if sign > 0:
            value += n//prod
        else:
            value -= n//prod

        for i in xrange(k, num):
            p_prod = primes[i]*prod
            if p_prod <= n:
                push((p_prod, i + 1, -1*sign))
            else:
                break

    return value


def mapes(n):
    root = int(n**(2.0/3.0))
    primes = rosemary.number_theory.sieves.primes(root)
    t = n**(0.33333333333333)

    c = bisect(primes, t)
    b = bisect(primes, n**(0.5))

    value = (b + c - 2)*(b - c + 1)//2

    for i in xrange(c, b):
        idx = bisect(primes, n//primes[i])
        value -= idx

    def P2(x, a):
        #total = 0
        #for i in xrange(a, len(primes)):
        #    for j in xrange(i, len(primes)):
        #        pp = primes[i - 1]*primes[j - 1]
        #        if pp <= x:
        #            total += 1

        #return total
        value = 0
        for i in xrange(a, len(primes)):
            p = primes[i - 1]

            if p*p > x:
                break

            pi = bisect(primes, x//p)
            pip = bisect(primes, p)
            value += (pi - pip + 1)

        return value

    def phi(x, a, cache={}):
        if (x, a) in cache:
            return cache[x, a]

        if a <= 5:
            if a == 1:
                value = x - x//2

            elif a == 2:
                value = x - x//2 - x//3 + x//6

            elif a == 3:
                value = x - x//2 - x//3 - x//5 + x//6 + x//10 + x//15 - x//30

            elif a == 4:
                value = (x - x//2 - x//3 - x//5 - x//7 +
                         x//6 + x//10 + x//14 + x//15 + x//21 + x//35 -
                         x//30 - x//42 - x//70 - x//105 +
                         x//210)
            elif a == 5:
                value = (x - x//2 - x//3 - x//5 - x//7 - x//11 +
                         x//6 + x//10 + x//14 + x//22 + x//15 + x//21 + x//33 + x//35 + x//55 + x//77 -
                         x//30 - x//42 - x//66 - x//70 - x//110 - x//154 - x//105 - x//165 - x//231 - x//385 +
                         x//210 + x//330 + x//462 + x//770 + x//1155 -
                         x//2310)

        elif a >= bisect(primes, x**(0.5)):
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

    value += phi(n, c)
    return value


def meissel_lehmer(n):
    root = int(n**(2.0/3.0))
    primes = rosemary.number_theory.sieves.primes(root)
    t = n**(0.33333333333333)

    c = bisect(primes, t)
    b = bisect(primes, n**(0.5))

    value = (b + c - 2)*(b - c + 1)//2

    for i in xrange(c, b):
        idx = bisect(primes, n//primes[i])
        value -= idx

    cache = {}

    def phi(x, a):
        if (x, a) in cache:
            return cache[x, a]

        if a <= 4:
            if a == 4:
                return ((x + 1)//2 - (x + 3)//6 - (x + 5)//10 + (x + 15)//30 - (x + 7)//14 + (x + 21)//42 + (x + 35)//70
                        - (x + 105)//210)
            elif a == 3:
                return (x + 1)//2 - (x + 3)//6 - (x + 5)//10 + (x + 15)//30
            elif a == 2:
                return (x + 1)//2 - (x + 3)//6
            elif a == 1:
                return (x + 1)//2
        else:
            val = phi(x, a - 1)
            if x >= primes[a - 1]:
                val -= phi(x//primes[a - 1], a - 1)

            cache[x, a] = val
            return val

    value += phi(n, c)
    return value


def lmo(x):
    root = int(x**(2.0/3.0))

    # print "sieving"

    primes = rosemary.number_theory.sieves.primes(root)
    t = x**(0.33333333333333)

    c = bisect(primes, t)
    b = bisect(primes, x**(0.5))

    value = (b + c - 2)*(b - c + 1)//2

    for i in xrange(c, b):
        idx = bisect(primes, x//primes[i])
        value -= idx

    # special = []
    # special_append = special.append
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

    # print "processing"

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


def mapes2(n):
    root = int(n**(2.0/3.0))
    primes = rosemary.number_theory.sieves.primes(root)
    t = n**(0.33333333333333)

    c = bisect(primes, t)
    b = bisect(primes, n**(0.5))

    value = (b + c - 2)*(b - c + 1)//2

    for i in xrange(c, b):
        idx = bisect(primes, n//primes[i])
        value -= idx

    cutoff = 5
    cache = {}
    mk = {}
    phi_mk = {}

    for k in xrange(2, cutoff + 1):
        mk[k], phi_mk[k], cache[k] = lookup[k]

    def phi(x, a):
        if x < primes[a - 1]:
            return 1

        elif a <= cutoff:
            if a == 1:
                return (x + 1)//2
            else:
                return phi_mk[a]*(x//mk[a]) + cache[a][x % mk[a]]

        else:
            if (x, a) in cache:
                return cache[x, a]

            elif a >= bisect(primes, x**(0.5)):
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

    value += phi(n, c)
    return value


if __name__ == "__main__":
    n = int(sys.argv[1])
    print lmo(n)
