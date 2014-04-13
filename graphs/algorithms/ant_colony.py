import math
import random
import numpy

def euc_2d(c1, c2):
    dist = math.sqrt((c1[0] - c2[0])**2 + (c1[1] - c2[1])**2)
    return dist

def cost(permutation, cities):
    dist = 0
    for (i, c1) in enumerate(permutation):
        c2 = permutation[0] if i == len(permutation) - 1 else permutation[i + 1]
        dist += euc_2d(cities[c1], cities[c2])
    return dist

def random_permutation(cities):
    new_cities = list(cities)
    random.shuffle(new_cities)
    return new_cities

def initialize_pheromone_matrix(num_cities, init_pher):
    A = numpy.ones((num_cities, num_cities))*init_pher
    return A

def calculate_choices(cities, last_city, exclude, pheromone, c_heur, c_hist):
    choices = []
    for (i, coord) in enumerate(cities):
        if i in exclude:
            continue

