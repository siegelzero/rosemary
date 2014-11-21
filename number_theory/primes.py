import rosemary.number_theory.sieves
from rosemary.number_theory.tables import lookup
import bisect

from collections import defaultdict

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


def meissel_lehmer(n):
    root = int(n**(2.0/3.0))
    primes = rosemary.number_theory.sieves.primes(root)

    c = sum(1 for p in primes if p*p*p <= n)
    b = sum(1 for p in primes if p*p <= n)

    pi = [0]*root
    prime_set = set(primes)
    count = 0

    for k in xrange(root):
        if k in prime_set:
            count += 1
        pi[k] = count

    value = (b + c - 2)*(b - c + 1)//2
    value -= sum(pi[n//primes[i]] for i in xrange(c, b))
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

    c = bisect.bisect(primes, t)
    b = bisect.bisect(primes, x**(0.5))

    value = (b + c - 2)*(b - c + 1)//2

    for i in xrange(c, b):
        idx = bisect.bisect(primes, x//primes[i])
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


def meissel_lehmer_lookup(x):
    root = int(x**(2.0/3.0))
    primes = rosemary.number_theory.sieves.primes(root)

    c = sum(1 for p in primes if p*p*p <= x)
    b = sum(1 for p in primes if p*p <= x)

    pi = [0]*root
    prime_set = set(primes)
    count = 0

    for k in xrange(root):
        if k in prime_set:
            count += 1
        pi[k] = count

    value = (b + c - 2)*(b - c + 1)//2
    value -= sum(pi[x//primes[i]] for i in xrange(c, b))

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
            return phi_mk[a]*(x//mk[a]) + cache[a][x % mk[a]]
        else:
            if a in cache:
                if x in cache[a]:
                    return cache[a][x]
            else:
                cache[a] = {}

            cache[a][x] = phi(x, a - 1) - phi(x//primes[a - 1], a - 1)
            return cache[a][x]

    value += phi(x, c)

    return value


if __name__ == "__main__":
    n = int(sys.argv[1])
    print lmo(n)
