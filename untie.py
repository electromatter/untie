#!/usr/bin/env python3

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

class Segment:
    def __init__(self, start, end, n):
        start, end = start % n, end % n
        self.start = start
        self.end = end
        self.n = n
        if start <= end:
            self.len = end - start + 1
        else:
            self.len = end - start + n + 1

    def __len__(self):
        return self.len

    def __call__(self, i):
        if i < -len(self) or i >= len(self):
            raise IndexError('index out of range')
        if i >= 0:
            return (self.start + i) % self.n
        else:
            return (self.end + i) % self.n

    def intersect(self, other):
        if self.n != other.n:
            raise ValueError('Fields do not match')
        dist = (other.start - self.start) % self.n
        return dist < self.len or self.n - dist < other.len

    def sub_seg(self, other):
        # This could seriously be reduced to a few conditionals
        segments = self.split()
        for cut_start, cut_end in other.split():
            done = []
            # Apply cuts
            for seg_start, seg_end in segments:
                if seg_start < cut_end < seg_end:
                    done.append((cut_end + 1, seg_end))
                if seg_start < cut_start < seg_end:
                    done.append((seg_start, cut_start - 1))
                if cut_end < seg_start or seg_end < cut_start:
                    done.append((seg_start, seg_end))
            segments = done

        # Wrap around
        done = []
        first, last = None, None
        for segment in segments:
            if segment[1] == self.n:
                last = segment
            elif segment[0] == 0:
                first = segment
            else:
                done.append(segment)

        if first and last:
            done.append((last[0], first[1]))
        elif first:
            done.append(first)
        elif last:
            done.append(last)
        segments = done

        return set(Segment(start, end, self.n) for start, end in segments)

    def split(self):
        if self.start < self.end:
            return [(self.start, self.end)]
        if self.end == 0:
            return [(self.start, self.n)]
        return [(self.start, self.n),
                (0, self.end)]

    def __repr__(self):
        return 'Segment(%i, %i, n=%i)' % (self.start, self.end, self.n)

def dplus(i, j, n):
    if i <= j:
        return j - i
    else:
        return j - i + n

class KnotSpinView(collections.abc.Sequence):
    def __init__(self, untie, segment):
        self.untie = untie
        self.segment = segment

    def __len__(self):
        return len(self.segment)

    def __getitem__(self, i):
        return self.untie.spin[self.segment(i)]

class Knot(collections.abc.MutableSequence):
    def __init__(self, untie, segment):
        self.untie = untie
        self.segment = segment
        self.step_count = 0
        self.spin = KnotSpinView(self.untie, self.segment)

        if len(self.segment) < 1:
            raise ValueError('Knots must be at least two elements')

    def __len__(self):
        return len(self.segment)

    def __getitem__(self, i):
        return self.untie.perm[self.segment(i)]

    def __setitem__(self, i, val):
        self.untie.perm[self.segment(i)] = val

    def __delitem__(self, i):
        raise NotImplementedError

    def insert(self, val, i=0):
        raise NotImplementedError

    def comparable(self, i, j):
        return self.untie.comparable(self.segment(i), self.segment(j))

    def swap(self, i, j):
        self.untie.swap(self.segment(i), self.segment(j))

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
        final = set()
        for seg in knots:
            final |= seg.sub_seg(self.segment)
        return final

    def __repr__(self):
        return '<Knot on %r step=%i data=%r spin=%r>' % (self.segment, self.step_count, list(self), list(self.spin))

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

    def add(self, pair):
        i, j = pair
        if i in self.nodes or j in self.nodes:
            raise ValueError('Adding pair to matching fails disjoint')
        self.pairs.add(pair)
        self.nodes.add(i)
        self.nodes.add(j)

    def discard(self, pair):
        i, j = pair
        self.pairs.discard(pair)
        self.nodes.discard(i)
        self.nodes.discard(j)

    def __repr__(self):
        return '<Matching pairs=%r>' % self.pairs

class Untie:
    def __init__(self, perm):
        self.perm = list(perm)
        if set(self.perm) != set(range(self.n)):
            raise ValueError('perm must be a permutation of [0, 1, ..., n-1]')
        self.spin = spin(self.perm)
        self.pairs = self.find_pairs()
        self.knots = []
        self.find_new_knots()
        self.matchings = []
        self.current_matching = Matching(self.n)

    @property
    def n(self):
        return len(self.perm)

    def find_knots(self):
        knots = []
        start = None
        for i in range(1, self.n):
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
        return set((i, j) for i, j in knots)

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

        if pair in self.current_matching:
            self.current_matching.remove(pair)
        else:
            self.current_matching.add(pair)

    @property
    def step_count(self):
        return len(self.matchings) + (len(self.current_matching) > 0)

    def finish_step(self):
        if self.current_matching:
            self.matchings.append(self.current_matching)
            self.current_matching = Matching(self.n)

    def find_new_knots(self):
        knots = [Segment(i, j, self.n) for i, j in self.find_knots()]
        for knot in self.knots:
            knots = knot.split(knots)
        knots = [knot for knot in knots if len(knot) > 1]
        for knot in knots:
            self.knots.append(Knot(self, knot))

    def step(self):
        if not self.knots:
            return False

        for knot in self.knots:
            knot.step()

        self.knots = [knot for knot in self.knots if not knot.is_sorted]

        # Find new knots
        self.find_new_knots()

        self.finish_step()
        return bool(self.pairs)

    def run(self):
        while self.step():
            pass
        return self.step_count

    def __repr__(self):
        fields = []
        fields.append('steps=%r' % self.step_count)
        fields.append('data=%r' % self.perm)
        fields.append('spin=%r' % self.spin)
        fields.append('pairs=%r' % self.pairs)
        fields.append('matchings=%r' % self.matchings)
        if self.current_matching:
            fields.append('current=%r' % self.current_matching)
        fields.append('knots=%r' % self.knots)
        return '<Untie %s>' % ' '.join(fields)

def main():
    import itertools

    N = 7

    for perm in itertools.permutations(range(N)):
        u = Untie(perm)
        if u.run() >= N:
            print(perm)

if __name__ == '__main__':
    main()
