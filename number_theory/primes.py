import rosemary.number_theory.sieves
from rosemary.number_theory.tables import lookup
from rosemary.data_structures import bit_sieve
from bisect import bisect

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


def lmo_bit(x):
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

    #block = bit_sieve.BITSieve(root)
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


if __name__ == "__main__":
    n = int(sys.argv[1])
    print lmo_bit(10**n)
    #print mapes(10**n)
