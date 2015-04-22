#!/usr/bin/pypy
import random


class TravelingSalesman(object):
    def __init__(self, points):
        self.points = points
        self.n = len(points)

        self.nearest_neighbors = []

        n = self.n
        k = min(n, 20)

        for i in xrange(n):
            neighbors = []
            xi, yi = points[i]
            for j in xrange(n):
                if i == j:
                    continue
                xj, yj = points[j]
                dist = ((xi - xj)**2 + (yi - yj)**2)**(0.5)
                neighbors.append((dist, j))

            neighbors.sort()
            self.nearest_neighbors.append(neighbors[0:k])

    def cost(self, path):
        points = self.points
        x0, y0 = points[path[0]]
        x1, y1 = points[path[-1]]
        total = ((x0 - x1)**2 + (y0 - y1)**2)**(0.5)
        for i in xrange(self.n - 1):
            x0, y0 = points[path[i]]
            x1, y1 = points[path[i + 1]]
            total += ((x0 - x1)**2 + (y0 - y1)**2)**(0.5)
        return total

    def two_opt(self, path):
        def rev(path, i, j):
            r = list(path)
            if i < j:
                r[i:j + 1] = r[i:j + 1][::-1]
            elif j < i:
                n = len(r)
                k = i
                idxs = []
                while True:
                    idxs.append(k)
                    k = (k + 1) % n
                    if k == j:
                        idxs.append(k)
                        break

                vals = [path[l] for l in idxs]
                rvals = vals[::-1]

                for (i, j) in enumerate(idxs):
                    r[j] = rvals[i]

            return r

        n = self.n
        points = self.points

        i = 0
        while i < n:
            j = (i + 1) % n
            xi = path[i]
            xj = path[j]
            (x0, y0) = points[xi]
            (x1, y1) = points[xj]
            dij = ((x0 - x1)**2 + (y0 - y1)**2)**(0.5)

            for (djk, xk) in self.nearest_neighbors[xj]:
                if djk >= dij:
                    break
                k = path.index(xk)
                l = (k - 1) % n
                xl = path[l]

                (x2, y2) = points[xl]
                (x3, y3) = points[xk]

                gain = dij
                gain += ((x2 - x3)**2 + (y2 - y3)**2)**(0.5)
                gain -= djk
                gain -= ((x0 - x2)**2 + (y0 - y2)**2)**(0.5)

                if gain > 0.000001:
                    path = rev(path, k, i)
                    # if i < k:
                    #     path[i:k + 1] = path[i:k + 1][::-1]
                    # else:
                    #     r_path = (path[k:] + path[:i + 1])[::-1]
                    #     path[k:] = r_path[k:]
                    #     path[:i + 1] = r_path[:i + 1]

                    i = -1
                    break

            i += 1

        return (self.cost(path), path)

    def select(self, popsize):
        print "Creating initial population"
        population = []

        for i in xrange(popsize):
            print "{}/{}".format(i + 1, popsize)
            path = range(self.n)
            random.shuffle(path)
            cost, opt = self.two_opt(path)
            population.append((cost, opt))

        return population

    def mgkrec(self, A, B):
        n = self.n
        h = random.randint(10, n//2)
        ret_vals = []

        for (P1, P2) in [(B, A), (A, B)]:
            j = random.randint(0, n - 1)
            C = [0]*n
            T = set()

            for i in xrange(h):
                C[i] = P1[(i + j) % n]
                T.add(C[i])

            i = h
            for p in P2:
                if p in T:
                    continue
                C[i] = p
                i += 1
                if i >= n:
                    break

            ret_vals.append(self.two_opt(C))

        return ret_vals

    def solve(self, popsize, cmax):
        population = self.select(popsize)

        for i in xrange(popsize):
            population.append(('inf', []))

        population.sort()

        best_cost, Xbest = population[0]

        for c in xrange(cmax):
            for i in xrange(popsize//2):
                p2i = population[2*i][1]
                p2i1 = population[2*i + 1][1]
                ((A_cost, A), (B_cost, B)) = self.mgkrec(p2i, p2i1)

                population[popsize + 2*i] = (A_cost, A)
                population[popsize + 2*i + 1] = (B_cost, B)

            population.sort()
            current_cost = population[0][0]

            if current_cost < best_cost:
                Xbest = list(population[1])
                best_cost = current_cost
                print c, best_cost

        return (best_cost, Xbest)


if __name__ == "__main__":
    berlin52 = [
        [565, 575], [25, 185], [345, 750], [945, 685], [845, 655], [880, 660], [25, 230], [525, 1000], [580, 1175],
        [650, 1130], [1605, 620], [1220, 580], [1465, 200], [1530, 5], [845, 680], [725, 370], [145, 665], [415, 635],
        [510, 875], [560, 365], [300, 465], [520, 585], [480, 415], [835, 625], [975, 580], [1215, 245], [1320, 315],
        [1250, 400], [660, 180], [410, 250], [420, 555], [575, 665], [1150, 1160], [700, 580], [685, 595], [685, 610],
        [770, 610], [795, 645], [720, 635], [760, 650], [475, 960], [95, 260], [875, 920], [700, 500], [555, 815],
        [830, 485], [1170, 65], [830, 610], [605, 625], [595, 360], [1340, 725], [1740, 245]
    ]

    t = TravelingSalesman(berlin52)
    print t.solve(16, 400)
