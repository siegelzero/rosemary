import math
import random
from collections import namedtuple
import rosemary.graphs.algorithms.tsp

def dist(c1, c2):
    d2 = (c1[0] - c2[0])**2 + (c1[1] - c2[1])**2
    return math.sqrt(d2)

def cost(permutation, cities):
    distance = 0
    for (i, c1) in enumerate(permutation):
        if i == len(permutation) - 1:
            c2 = permutation[0]
        else:
            c2 = permutation[i + 1]
        distance += dist(cities[c1], cities[c2])
    return distance

def random_two_opt(current):
    path = current.vector
    num_points = len(path)
    c1 = random.randint(0, num_points - 1)
    c2 = random.randint(0, num_points - 1)
    exclude = set([c1])

    if c1 == 0:
        exclude.add(num_points - 1)
    else:
        exclude.add(c1 - 1)

    if c1 == num_points - 1:
        exclude.add(0)
    else:
        exclude.add(c1 + 1)

    if c2 < c1:
        c1, c2 = c2, c1

    perm = list(path)
    perm[c1:c2] = perm[c1:c2][::-1]
    return perm

def create_neighbor(current, cities):
    candidate_vector = random_two_opt(current)
    candidate_cost = cost(candidate_vector, cities)
    candidate = namedtuple('Candidate', ('cost', 'vector'))
    candidate.cost = candidate_cost

    candidate.vector = candidate_vector

    return candidate

def should_accept(candidate, current, temp):
    if candidate.cost < current.cost:
        return True
    else:
        return math.exp((current.cost - candidate.cost) / temp) > random.random()

def search(cities, max_iter=None, max_temp=None, temp_change=None):
    if max_iter is None:
        max_iter = 200000
    if max_temp is None:
        max_temp = 100000.0
    if temp_change is None:
        temp_change = 0.9999

    current = namedtuple('Candidate', ('cost', 'vector'))
    current.vector = range(len(cities))
    random.shuffle(current.vector)
    current.cost = cost(current.vector, cities)
    #c, v = rosemary.graphs.algorithms.tsp.steepest_ascent_two_opt(cities)
    #current.vector = v
    #current.cost = c
    #print c

    temp = max_temp
    best = current

    for j in xrange(10):
        print j
        for i in xrange(max_iter):
            candidate = create_neighbor(current, cities)
            temp *= temp_change

            if should_accept(candidate, current, temp):
                current = candidate

            if candidate.cost < best.cost:
                best = candidate
                print best.cost

        current = best
        temp = max_temp

    return best

if __name__ == "__main__":
    berlin52 = [[565,575],[25,185],[345,750],[945,685],[845,655],[880,660],[25,230],[525,1000],[580,1175],[650,1130],[1605,620],[1220,580],[1465,200],[1530,5],[845,680],[725,370],[145,665],[415,635],[510,875],[560,365],[300,465],[520,585],[480,415],[835,625],[975,580],[1215,245],[1320,315],[1250,400],[660,180],[410,250],[420,555],[575,665],[1150,1160],[700,580],[685,595],[685,610],[770,610],[795,645],[720,635],[760,650],[475,960],[95,260],[875,920],[700,500],[555,815],[830,485],[1170,65],[830,610],[605,625],[595,360],[1340,725],[1740,245]]
    max_iterations = 150000
    max_temp = 100000.0
    temp_change = 0.9999

    best = search(berlin52, max_iterations, max_temp, temp_change)

    print "Done. Best Solution: {}, {}".format(best.cost, best.vector)

