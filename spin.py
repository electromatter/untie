#!/usr/bin/env python3

def displacement(perm):
    perm = list(perm)
    if set(perm) != set(range(len(perm))):
        raise ValueError('must be a permutation of [0, 1, 2, ..., n-1]')
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

def main():
    while True:
        try:
            perm = eval(input())
            print(perm, '=>', spin(perm))
        except Exception:
            print()
        except (KeyboardInterrupt, EOFError):
            break

if __name__ == '__main__':
    main()
