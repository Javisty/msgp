'''
A set of utility objects for the MSGP algorithms.
'''
import numpy as np
from os import listdir
from os.path import isfile, join

def get_gradual_items(data):
    '''Returns wether the indices increases, decreases or fluctuate
    through first dimension of data.
    data should be a 2D numpy array.
    Output is a numpy array of size the second dimension of data.
    Output[i] is 1/0/-1 if column i has its values increasing/
    fluctuating/decreasing through data.'''

    increases = np.diff(data, axis=0) >= 0
    return np.array([1 if increases[:,i].all()
                     else -1 if (~increases[:,i]).all()
                     else 0 for i in range(increases.shape[1])], dtype=np.int8)

def remove_non_frequent(gradualities, minSup):
    '''gradualities[s, m, i] returns the graduality (1, -1 or 0) of attribute i
    between stages s and s+1 in cycle m.
    Changes inplace to put to 0 non-frequent attributes for each starting stage s.'''
    M = gradualities.shape[1]

    # For each stage and each attribute, find number of increase over cycles
    up_occurences = np.count_nonzero(gradualities == 1, axis=1)

    up_non_frequent = np.swapaxes(np.tile(up_occurences < minSup, (M, 1, 1)), 0, 1)
    down_non_frequent = np.swapaxes(np.tile(M-up_occurences < minSup, (M, 1, 1)), 0, 1)

    gradualities[up_non_frequent & (gradualities == 1)] = 0
    gradualities[down_non_frequent & (gradualities == -1)] = 0

def print_results_seasons(results, k, I):
    print("Season -> gradual patterns")
    for start, length in results.items():
        for l, patterns in length.items():
            print("From", start, "to", (start+l)%k, "->", end=' ')
            for pattern in patterns:
                for item in pattern:
                    if item < I:  # increasing item
                        print(str(item) + str('+'), end='')
                    else:  # decreasing item
                        print(str(item-I) + str('-'), end='')
                    print(', ', end='')
            print()

def print_results_patterns(results, k, I):
    print("Pattern -> seasons")
    for pattern, seasons in results.items():
        print(*[str(i)+"+" if i<I else str(i-I)+"-" for i in pattern],
              sep=',', end=" -> ")
        for start, length in seasons:
            print(start, "-", (start+length)%k, sep='', end=', ')
        print()

def check_pattern(season, pattern, data, minSup):
    '''Check if pattern is respected through period over data, w.r.t. minSup.
    pattern (set): contains the gradual items encoded from 0 to 2*I-1
    season: tuple, first element is the starting index of the season,
    second element is the length 0<l<k+1 of the season.'''
    M, k, I = data.shape
    s, l = season
    increasing = [ i for i in list(pattern) if i<I]
    decreasing = [ i-I for i in list(pattern) if i>=I]

    if s+l < k:  # not cross-cycle season
        respected_incr = np.all((np.diff(data[:,s:s+l+1,increasing], axis=1) >= 0), axis=(1,2))
        respected_decr = np.all((np.diff(data[:,s:s+l+1,decreasing], axis=1) < 0), axis=(1,2))

    else:  # remove last cycle
        season_data = data.reshape((M*k, I))[s:-(k-s),:].reshape((M-1,k,I))
        respected_incr = np.all(np.diff(season_data[:,:l+1,increasing], axis=1) >= 0, axis=(1,2))
        respected_decr = np.all(np.diff(season_data[:,:l+1,decreasing], axis=1) < 0, axis=(1,2))

    if not increasing:  # no increase gradual item in pattern
        return np.count_nonzero(respected_decr) >= minSup
    elif not decreasing:
        return np.count_nonzero(respected_incr) >= minSup
    else:
        return np.count_nonzero(respected_incr & respected_decr) >= minSup

def check_results_seasons(data, results, minSup):
    failed = set()
    for start, length in results.items():
        for l, patterns in length.items():
            for pattern in patterns:
                if not check_pattern((start, l), pattern, data, minSup):
                    failed.add((start, l, pattern))
    print("Numbers of fail:", len(failed))
    return failed

def check_results_patterns(data, results, minSup):
    failed = set()
    for pattern, seasons in results.items():
        for season in seasons:
            start, l = season
            if not check_pattern((start, l), pattern, data, minSup):
                failed.add((start, l, pattern))
    print("Numbers of fail:", len(failed))
    return failed

def check_testcases(dirpath, algo, minSup):
    failed = set()
    for testcase in listdir(dirpath):
        if not "50.npy" in testcase:
            print("Testing", testcase)
            data = np.load(join(dirpath, testcase))
            results = algo(data, minSup)
            failure_testcase = check_results_periods(data, results, minSup)
            if failure_testcase:
                print("Testcase", testcase, "failed")
                failed.add(failure_testcase)
    return failed

def compare_results_seasons(res1, res2):
    '''Returns the differences between the gradual patterns
    stored in res1 and res2, in the season format.'''
    diff = set()
    for start, length in res1.items():
        for l, patterns in length.items():
            if (res2.get(start)) and not (patterns == res2[start].get(l)):
                diff.add((start, l))
    return diff

def compare_results_patterns(res1, res2):
    '''Returns the differences between the gradual patterns
    stored in res1 and res2, in the pattern format.'''
    diff = set()
    for pattern, seasons in res1.items():
        seasons2 = res2.get(pattern)
        if not(seasons2) and seasons:
            diff.add(pattern)
        elif seasons != seasons2:
            diff.add(pattern)
    return diff

def count_patterns_seasons(patterns):
    '''Return the number of patterns contained in patterns,
    with the result format of MSGP_season'''
    return sum([ sum([ len(k) for k in start.values() ])
                  for start in patterns.values()])

def count_patterns_patterns(patterns):
    '''Return the number of patterns contained in patterns,
    with the result format of MSGP_patterns'''
    return sum([ len(seasons) for seasons in patterns.values()])

def cover(i, data):
    '''Returns the list of timestamps such that i is respected on the
    following transition.
    i is a gradual item, <I if increasing, >=I otherwise.'''
    _, I = data.shape

    if i < I:  # increase
        return np.where(np.diff(data[:,i]) >= 0)[0]
    else:  # decrease
        return np.where(np.diff(data[:,i-I]) < 0)[0]
