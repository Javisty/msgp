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

def print_results_seasons(results, k):
    print("Season -> gradual patterns")
    for start, length in results.items():
        for l, patterns in length.items():
            print(''.join(map(str, range(start, start+length))), ":", end=' ')
            for pattern in patterns:
                print('+'.join(map(str, pattern[0])) + '-'.join(map(str, pattern[1])),
                      end=', ')
            print()

def check_pattern(season, pattern, data, minSup):
    '''Check if pattern is respected through period over data, w.r.t. minSup.
    pattern: tuple, first set contains indices of increasing attributes,
    second set contains indices of decreasing attributes.
    season: tuple, first element is the starting index of the season,
    second element is the length 0<l<k+1 of the season.'''
    i, l = season
    increasing, decreasing = list(pattern[0]), list(pattern[1])

    if i+l < data.shape[1]:
        count_increasing = np.sum((np.diff(data[:,i:i+l,increasing], axis=1) >= 0).all(axis=(1,2)))
        count_decreasing = np.sum((np.diff(data[:,i:i+l,decreasing], axis=1) < 0).all(axis=(1,2)))
    else:
        count_increasing = np.sum((np.diff(data[:-1,i:i+l,increasing], axis=1) >= 0).all(axis=(1,2)))
        count_decreasing = np.sum((np.diff(data[:-1,i:i+l,decreasing], axis=1) < 0).all(axis=(1,2)))

    return (count_increasing + count_decreasing) >= minSup

def check_results_seasons(data, results, minSup):
    failed = set()
    for start, length in results.items():
        for l, patterns in length.items():
            for pattern in patterns:
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

def count_patterns_seasons(patterns):
    '''Return the number of patterns contained in patterns,
    with the result format of MSGP_season'''
    return sum([ sum([ len(k) for k in start.values() ])
                  for start in patterns.values()])

def cover(i, data):
    '''Returns the list of timestamps such that i is respected on the
    following.
    i is a gradual item, encoded as (attribute, bool) where bool is True
    when the increase is considered'''
    M, _ = data.shape

    return np.where((np.diff(data[:,i[0]]) >= 0) == i[1])[0]
