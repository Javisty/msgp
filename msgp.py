'''
Implementation of the second proposal of MSGP
'''
import numpy as np
from itertools import chain, combinations
import utils


def MFI(data, minSup):
    '''Implementation of the Eclat algorithm, for gradual signature.
    data is a gradual signature, of shape (M, I), with -1/0/1 values'''
    _, I = data.shape

    # Find frequent gradual items
    increasing_items = np.where(np.count_nonzero(data == 1, axis=0) >= minSup)[0]
    decreasing_items = np.where(np.count_nonzero(data == -1, axis=0) >= minSup)[0]

    # gradual itemsets are indexed from 0 to (2*I-1)
    itemsets = [ (frozenset({i}), frozenset(np.where(data[:,i]==1)[0])) for i in increasing_items ]\
        + [ (frozenset({I+i}), frozenset(np.where(data[:,i]==-1)[0])) for i in decreasing_items ]

    results = set([ itemset for itemset, _ in itemsets])
    l = 2
    while itemsets and l <= I:
        new_itemsets = list()
        for idx, (itemset, tids) in enumerate(itemsets):
            for idx2 in range(idx+1, len(itemsets)):
                itemset2, tids2 = itemsets[idx2]
                if sorted(list(itemset))[:l-2] == sorted(list(itemset2))[:l-2]:
                    intersect = tids.intersection(tids2)
                    if len(intersect) >= minSup:
                        new_itemset = itemset.union(itemset2)
                        new_itemsets.append((new_itemset, intersect))
                        results.add(new_itemset)
        itemsets = new_itemsets
        l += 1
    return results

def intersect(data1, data2):
    '''Returns the intersected transaction database.'''
    res = data1.copy()
    res[data1 != data2] = 0
    return res

def MFCCP(Gamma, k, minSup, verbose):
    # In case of verbose
    verboseprint = print if verbose else lambda *a, **k: None

    res = dict([(i, dict()) for i in range(k)])
    M, I = Gamma.shape[1:]
    l = 1
    to_visit = range(k)

    utils.remove_non_frequent(Gamma, minSup)
    while to_visit and l <= k:
        Lambda = list()
        for i in to_visit:
            patterns = MFI(Gamma[i], minSup)
            if patterns:
                Lambda.append(i)
                res[i][l] = patterns.copy()
            if i == k-1:  # handle inter-season difficulty
                Gamma[i] = intersect(Gamma[i], np.vstack((Gamma[0][1:,:], np.zeros((1,I)))))
            else:
                Gamma[i] = intersect(Gamma[i], Gamma[i+1])
        verboseprint("Patterns for periods of size", l, "found.", len(to_visit),
              "candidates to visit for size ", l+1, end='\r')
        l +=1
        to_visit = Lambda[:]

    return res

def MSGP_seasons(data, minSup, verbose=True):
    '''
    Mining Seasonal Gradual Pattern. Returns all the frequent seasonal
    gradual patterns in data.
    Attributes and timestamps will be considered as indices, the decoding
    is left to the user.

    Parameters:
    -----------
    data (numpy.array): an array of shape (M, k, I) with numerical values.
    First dimension stands for the seasons, second for the observations,
    last one for the attribute values. The database should be complete.

    minSup (integer): the minimum number of seasons where a pattern
    should occur to be frequent.

    Output:
    -------
    delta: a dictionary of all frequent seasonal gradual patterns
    '''
    # In case of verbose
    verboseprint = print if verbose else lambda *a, **k: None

    # Number of seasons, season length and number of attributes
    M, k, I = data.shape

    Gamma = np.zeros((k, M, I), dtype=np.int8)
    verboseprint("Initialization...", end='\r')
    for stage in range(k):
        for season in range(M):
            if stage != k-1:
                Gamma[stage, season] = utils.get_gradual_items(data[season, [stage, (stage+1)%k]])
            else:  # special treatment for last stage
                if season != M-1:  # use first stage of next season
                    Gamma[stage, season] = utils.get_gradual_items(np.vstack((data[season, stage, :], data[season+1, 0, :])))
                else:  # for last stage of last season, no data
                    Gamma[stage, season] = np.zeros(I)

    verboseprint("Starting BFS", end='\r')
    return MFCCP(Gamma, k, minSup, verbose)

##################### Second algorithm #####################
# UNCOMPLETE

def CFP(w, N, m):
    pass

def MFSP(Gamma, S, minSup):
    I = len(Gamma.keys())//2
    k = 1

    candidates = set([(frozenset([i]), frozenset()) for i in range(I)])
    candidates = candidates.update([(frozenset(), frozenset([i])) for i in range(I)])

    results = dict()
    while candidates and k <= I:
        psi = set()
        for p in candidates:
            modulo = (np.where(Gamma[p])[0])%S
            periods = CFS(modulo, S, minSup)
            if periods:
                results[p] = periods



def MSGP_patterns(data, k, minSup):
    '''
    Mining Seasonal Gradual Pattern. Returns all the frequent seasonal
    gradual patterns in data.
    Attributes and timestamps will be considered as indices, the decoding
    is left to the user.

    Parameters:
    -----------
    data (numpy.array): an array of shape (N=M*k, I) with numerical values.
    First dimension stands for the timestamps (k by cycle), second one for
    the attribute values. The database should be complete.

    k (int): the cycle-length of the data.

    minSup (integer): the minimum number of seasons where a pattern
    should occur to be frequent.

    Output:
    -------
    delta (list): list of the frequent seasonal patterns.
    '''
    M, I = data.shape

    Gamma = dict()
    for i in range(I):
        grad_incr = (frozenset([i]), frozenset())
        grad_decr = (frozenset(), frozenset([i]))

        Gamma[grad_incr] = utils.cover((i, True), data)
        Gamma[grad_decr] = utils.cover((i, False), data)

    SI = MFSP(Gamma, S, minSup)

    pass

##################### Brute Force #####################
# Test every possible seasonal gradual pattern
# Used to validate the results of the two previous algorithms

def brute_force(data, k, minSup, format='season'):
    '''
    Returns all the frequent seasonal gradual patterns in data.
    Brute Force proceeds by scanning the database for each seasonal
    gradual pattern. Very computational-heavy, however we are sure
    the results are good.
    Attributes and timestamps will be considered as indices, the decoding
    is left to the user.

    Parameters:
    -----------
    data (numpy.array): an array of shape (N=M*k, I) with numerical values.
    First dimension stands for the timestamps (k by cycle), second one for
    the attribute values. The database should be complete.

    k (int): the cycle-length of the data.

    minSup (integer): the minimum number of seasons where a pattern
    should occur to be frequent.

    format ('season' or 'pattern'): choose the format of the output structure.
    'season' is same as the output of MSGP_seasons, 'pattern' is like
    MSGP_patterns.

    Output:
    -------
    delta (list): list of the frequent seasonal patterns.
    '''
    N, I = data.shape
    M = N//k

    if format == 'season':
        results = dict()
        for start in range(k): # enumerate all seasons
            results[start] = dict()
            for length in range(1, k+1):
                results[start][length] = set()

                # enumerate all gradual patterns
                grad_patterns = chain.from_iterable(combinations(range(2*I), l)
                                                    for l in range(1, 2*I+1))

                # Check each pattern over the current season
                for pattern in grad_patterns:
                    if utils.check_pattern((start, length), set(pattern),
                                           data.reshape((M,k,I)), minSup):
                        results[start][length].add(frozenset(pattern))

    return results
