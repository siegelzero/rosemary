import rosemary.number_theory.sieves


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
    print "Sieving..."
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

    print "Recursing..."

    def phi(x, a):
        if (x, a) in cache:
            return cache[x, a]

        if x == 0:
            val = 0
        elif a == 1:
            val = (x + 1)//2
        else:
            val = phi(x, a - 1) - phi(x//primes[a - 1], a - 1)

        cache[x, a] = val
        return val

    value += phi(n, c)
    print "Done. Cache size: {}".format(len(cache))
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

        if x == 0:
            val = 0
        elif a == 1:
            val = (x + 1)//2
        else:
            val = phi(x, a - 1) - phi(x//primes[a - 1], a - 1)

        cache[x, a] = val
        return val

    phi_table = {}
    cutoff = 7
    for k in xrange(3, cutoff):
        print k
        mk = 1
        for i in xrange(k):
            mk *= primes[i]

        mk2 = mk//2
        block = [1]*(mk2 + 1)
        block[0] = 0

        p = 2
        count = 1
        while count <= k:
            for i in xrange(p, mk2, p):
                block[i] = 0

            p += 1
            while block[p] == 0:
                p += 1
            count += 1

        table = []
        count = 0
        for i in xrange(1, mk2):
            if block[i] == 1:
                count += 1
                table.append((i, count))

        phi_table[k] = table

    return phi_table


def phi_table(k):
    primes = rosemary.number_theory.sieves.primes(100)
    mk = 1
    for i in xrange(k):
        mk *= primes[i]

    mk2 = mk//2
    block = [1]*mk2
    block[0] = 0

    p = 2
    count = 1
    while count <= k:
        for i in xrange(p, mk2, p):
            block[i] = 0

        p += 1
        while block[p] == 0:
            p += 1
        count += 1

    table = []
    count = 0
    for i in xrange(1, mk2):
        if block[i] == 1:
            count += 1
            table.append((i, count))

    return table
