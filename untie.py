import collections.abc

def displacement(perm):
    return [p_i - i for i, p_i in enumerate(perm)]

def spin(perm):
    disp = displacement(perm)
    while True:
        maximum, minimum = max(disp), min(disp)
        if maximum - minimum <= len(disp):
            break
        xi, ni = disp.index(maximum), disp.index(minimum)
        disp[xi] -= len(disp)
        disp[ni] += len(disp)
    return disp

def dplus(i, j, n):
    if i <= j:
        return j - i
    else:
        return j - i + n

class Knot(collections.abc.MutableSequence):
    def __init__(self, untie, i, j):
        self.untie = untie
        self.i, self.j = i, j
        self.step_count = 0

    def _real_index(self, i):
        if i < -len(self) or i >= len(self):
            raise IndexError('index out of range')
        return (i + self.i) % self.untie.n

    def __len__(self):
        return self.untie.dplus(self.i, self.j) + 1

    def __getitem__(self, i):
        return self.untie.perm[self._real_index(i)]

    def __setitem__(self, i, val):
        self.untie.perm[self._real_index(i)] = val

    def __delitem__(self, i):
        raise NotImplementedError

    def insert(self, val, i=0):
        raise NotImplementedError

    def comparable(self, i, j):
        return self.untie.comparable(self._real_index(i), self._real_index(j))

    def swap(self, i, j):
        self.untie.swap(self._real_index(i), self._real_index(j))

    @property
    def is_sorted(self):
        return self.step_count >= len(self)

    def step(self):
        if self.is_sorted:
            return False
        i = self.step_count % 2
        while i < len(self) - 1:
            if self.comparable(i, i + 1):
                self.swap(i, i + 1)
            i += 2
        self.step_count += 1
        return not self.is_sorted

    def split(self, knots):
        return set()

class Matching(collections.abc.MutableSet):
    def __init__(self, n):
        self.n = n
        self.pairs = set()
        self.nodes = set()

    def __len__(self):
        return len(self.pairs)

    def __contains__(self, pair):
        return pair in self.pairs

    def __iter__(self):
        return iter(self.pairs)

    def _unpack(self, pair):
        i, j = pair
        if abs(dplus(i, j, self.n)) != 1:
            raise ValueError('pair is not an adjacent pair')
        return i % self.n, j % self.n

    def add(self, pair):
        i, j = self._unpack(pair)
        if i in self.nodes or j in self.nodes:
            raise ValueError('Adding pair to matching fails disjoint')
        self.pairs.add(pair)
        self.nodes.add(i)
        self.nodes.add(j)

    def discard(self, pair):
        i, j = self._unpack(pair)
        self.pairs.discard(pair)
        self.nodes.discard(i)
        self.nodes.discard(j)

class Untie:
    def __init__(self, perm):
        self.perm = list(perm)
        if set(self.perm) != set(range(self.n)):
            raise ValueError('perm must be a permutation of [0, 1, ..., n-1]')
        self.spin = spin(self.perm)
        self.pairs = self.find_pairs()
        self.knots = [Knot(self, i, j) for i, j in self.find_knots()]
        self.matchings = []
        self.current_matching = Matching(self.n)

    @property
    def n(self):
        return len(self.perm)

    def find_knots(self):
        knots = []
        start = None
        for i in range(self.n):
            j = (i + 1) % self.n
            if self.spin[i] > self.spin[j] and start is None:
                # Keep track of decreasing spin
                start = i
            elif self.spin[i] <= self.spin[j] and start is not None:
                # Non decreasing spin
                knots.append((start, i))
                start = None
        else:
            if start is not None and knots and knots[0][0] == 0:
                # Merge knots that wrap around
                knots[0] = (start, knots[0][1])
        return set(knots)

    def dplus(self, i, j):
        return dplus(i, j, self.n)

    def comparable(self, i, j):
        return self.spin[i] - self.spin[j] > self.dplus(i, j)

    def find_pairs(self):
        pairs = set()
        for i in range(self.n):
            for j in range(self.n):
                if self.comparable(i, j):
                    pairs.add((self.perm[i], self.perm[j]))
        return pairs

    def swap(self, i, j):
        if abs(self.dplus(i, j)) > 1:
            raise ValueError('Tried to swap to a non-adjacent pair')
        elif self.dplus(i, j) == -1:
            i, j = j, i

        pair = (self.perm[i], self.perm[j])

        if pair in self.pairs:
            self.pairs.remove(pair)
        else:
            self.pairs.add(pair)

        self.perm[i], self.perm[j] = self.perm[j], self.perm[i]
        self.spin[i], self.spin[j] = self.spin[j] + 1, self.spin[i] - 1

        if (i, j) in self.current_matching:
            self.current_matching.remove((i, j))
        else:
            self.current_matching.add((i, j))

    @property
    def step_count(self):
        return len(self.matchings) + (len(self.current_matching) > 0)

    def finish_step(self):
        if self.current_matching:
            self.matchings.append(self.current_matching)
            self.current_matching = Matching(self.n)

    def find_new_knots(self):
        knots = self.find_knots()
        for knot in self.knots:
            knots = knot.split(knots)
        for i, j in knots:
            self.knots.append(Knot(self, i, j))

    def step(self):
        for knot in self.knots:
            knot.step()

        self.knots = [knot for knot in self.knots if not knot.is_sorted]

        # Find new knots
        self.find_new_knots()

        self.finish_step()
        return bool(self.knots)
