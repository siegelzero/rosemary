import rosemary.number_theory.sieves

from collections import defaultdict


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

        if a == 1:
            val = (x + 1)//2
        else:
            val = phi(x, a - 1) - phi(x//primes[a - 1], a - 1)

        cache[x, a] = val
        return val

    value += phi(n, c)
    return value


def lmo(x):
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

    special = []

    stack = [(1, c, 1)]
    push = stack.append
    pop = stack.pop

    while stack:
        (n, a, sign) = pop()

        if a == 1 and n*n*n <= x:
            if sign > 0:
                value += (x//n + 1)//2
            else:
                value -= (x//n + 1)//2
        elif n*n*n > x:
            special.append((n, a, sign))
        else:
            push((n, a - 1, sign))
            push((n*primes[a - 1], a - 1, -sign))

    grouped = defaultdict(list)
    while special:
        (n, a, sign) = special.pop()
        grouped[a].append((x//n, sign))

    for a in grouped:
        grouped[a].sort(reverse=True)

    groups = grouped.items()
    groups.sort(reverse=True)

    block = [1]*(root + 1)
    block[0] = 0

    while groups:
        (a, values) = groups.pop()
        p = primes[a - 1]

        block[p::p] = [0]*(root//p)
        count = 0
        v, sign = values.pop()

        for i in xrange(1, root + 1):
            if block[i]:
                count += 1

            if i == v:
                if sign > 0:
                    value += count
                else:
                    value -= count

                if values:
                    v, sign = values.pop()
                else:
                    break

    return value


def lmo2(x):
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

    special = []

    stack = [(1, c, 1)]
    push = stack.append
    pop = stack.pop

    while stack:
        (n, a, sign) = pop()

        if a == 1 and n*n*n <= x:
            if sign > 0:
                value += (x//n + 1)//2
            else:
                value -= (x//n + 1)//2
        elif n*n*n > x:
            special.append((n, a, sign))
        else:
            push((n, a - 1, sign))
            push((n*primes[a - 1], a - 1, -sign))

    cache = {}

    def phi(n, a):
        if a <= 4:
            if a == 4:
                return ((n + 1)//2 - (n + 3)//6 - (n + 5)//10 + (n + 15)//30 - (n + 7)//14 + (n + 21)//42 + (n + 35)//70
                        - (n + 105)//210)
            elif a == 3:
                return (n + 1)//2 - (n + 3)//6 - (n + 5)//10 + (n + 15)//30
            elif a == 2:
                return (n + 1)//2 - (n + 3)//6
            elif a == 1:
                return (n + 1)//2
        else:
            if (n, a) in cache:
                return cache[n, a]

            val = phi(n, a - 1) - phi(n//primes[a - 1], a - 1)

            cache[n, a] = val
            return val

    special.sort(key=lambda (n, a, s): a)

    for (n, a, s) in special:
        if s > 0:
            value += phi(x//n, a)
        else:
            value -= phi(x//n, a)

    return value
