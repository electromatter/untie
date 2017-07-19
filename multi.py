#!/usr/bin/env python3

import itertools
import multiprocessing

import untie

N = 8
THREADS = 4

def do_permutation(perm):
    u = untie.Untie(perm)
    u.run()
    return perm, u

def handle_result(result):
    perm, u = result
    if u.step_count >= N - 1:
        print(perm, u.step_count)

def main():
    with multiprocessing.Pool(THREADS) as pool:
        for result in pool.imap_unordered(do_permutation, itertools.permutations(range(N))): 
            handle_result(result)

if __name__ == '__main__':
    main()
