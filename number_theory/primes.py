import rosemary.number_theory.sieves

from rosemary.number_theory.tables import lookup

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


def meissel(n):
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

    value = (b + c - 2)*(b - c + 1)//2 - sum(pi[n//primes[i]] for i in xrange(c, b))

    stack = [(1, 0, 1)]
    push = stack.append
    pop = stack.pop

    while stack:
        (prod, k, sign) = pop()

        if sign > 0:
            value += n//prod
        else:
            value -= n//prod

        for i in xrange(k, c):
            if primes[i]*prod <= n:
                push((primes[i]*prod, i + 1, -sign))
            else:
                break

    return value


def meissel2(n):
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

    cutoff = 5
    cache = {}
    mk = {}
    phi_mk = {}

    for k in xrange(2, cutoff + 1):
        mk[k], phi_mk[k], cache[k] = lookup[k]

    def phi(x, a):
        # if a == 1:
        #     return (x + 1)//2
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

    value += phi(n, c)

    return value


def meissel_lehmer_nonrec(n):
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

    cutoff = 5
    cache = {}
    mk = {}
    phi_mk = {}

    for k in xrange(2, cutoff + 1):
        mk[k], phi_mk[k], cache[k] = lookup[k]

    stack = [(n, c, 1)]
    push = stack.append
    pop = stack.pop

    while stack:
        (x, a, sign) = pop()

        if a == 1:
            if sign > 0:
                value += (x + 1)//2
            else:
                value -= (x + 1)//2
        elif a <= cutoff:
            if sign > 0:
                value += phi_mk[a]*(x//mk[a]) + cache[a][x % mk[a]]
            else:
                value -= phi_mk[a]*(x//mk[a]) + cache[a][x % mk[a]]
        else:
            push((x, a - 1, sign))

            if primes[a - 1] <= x:
                push((x//primes[a - 1], a - 1, -sign))

    return value


def meissel3(n):
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

        if x == 0:
            return 0
        elif a == 1:
            return (x + 1)//2
        else:
            cache[x, a] = phi(x, a - 1) - phi(x//primes[a - 1], a - 1)
            return cache[x, a]

    value += phi(n, c)

    return value


def phi_table(k):
    primes = rosemary.number_theory.sieves.primes(100)
    mk = 1
    for i in xrange(k):
        mk *= primes[i]

    block = [1]*mk
    block[0] = 0

    p = 2
    count = 1
    while count <= k:
        for i in xrange(p, mk, p):
            block[i] = 0

        p += 1
        while block[p] == 0:
            p += 1
        count += 1

    table = [0]*mk
    count = 0

    for i in xrange(1, mk):
        if block[i] == 1:
            count += 1
        table[i] = count

    return table


if __name__ == "__main__":
    n = int(sys.argv[1])
    sys.setrecursionlimit(1500)
    print meissel2(n)
