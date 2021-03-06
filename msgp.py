'''
Implementation of the second proposal of MSGP
'''
import numpy as np
from itertools import chain, combinations
from copy import deepcopy
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

def MFCI(data, minSup):
    '''Implementation of the Eclat-Close algorithm , for gradual signature.
    data is a gradual signature, of shape (M, I), with -1/0/1 values'''
    _, I = data.shape

    # Find frequent gradual items
    increasing_items = np.where(np.count_nonzero(data == 1, axis=0) >= minSup)[0]
    decreasing_items = np.where(np.count_nonzero(data == -1, axis=0) >= minSup)[0]

    # gradual itemsets are indexed from 0 to (2*I-1)
    itemsets = [ (frozenset({i}), frozenset(np.where(data[:,i]==1)[0])) for i in increasing_items ]\
        + [ (frozenset({I+i}), frozenset(np.where(data[:,i]==-1)[0])) for i in decreasing_items ]

    # mapping t-IDset to closure set
    hash_table = {tidset: item for item, tidset in itemsets}

    l = 2
    while itemsets and l <= I:
        new_itemsets = list()
        for idx, (itemset, tids) in enumerate(itemsets):
            for idx2 in range(idx+1, len(itemsets)):
                itemset2, tids2 = itemsets[idx2]
                if sorted(list(itemset))[:l-2] == sorted(list(itemset2))[:l-2]:
                    intersect_tid = tids.intersection(tids2)
                    union_items = itemset.union(itemset2)
                    if len(intersect_tid) >= minSup:
                        if intersect_tid in hash_table:
                            hash_table[intersect_tid] = hash_table[intersect_tid].union(union_items)
                        else:
                            hash_table[intersect_tid] = union_items

                        new_itemsets.append((union_items, intersect_tid))
        itemsets = new_itemsets
        l += 1
    return set(hash_table.values())

def intersect(data1, data2):
    '''Returns the intersected transaction database.'''
    res = data1.copy()
    res[data1 != data2] = 0
    return res

def MFCCP(Gamma, k, minSup, closed, verbose):
    # In case of verbose
    verboseprint = print if verbose else lambda *a, **k: None
    # Frequent or Closed Frequent itemsets
    pattern_miner = MFCI if closed else MFI

    res = dict([(i, dict()) for i in range(k)])
    M, I = Gamma.shape[1:]
    l = 1
    to_visit = range(k)

    utils.remove_non_frequent(Gamma, minSup)
    while to_visit and l <= k:
        psi = deepcopy(Gamma)
        Lambda = list()
        for i in to_visit:
            patterns = pattern_miner(Gamma[i], minSup)
            if patterns:
                Lambda.append(i)
                res[i][l] = patterns.copy()
            if i == k-1:  # handle inter-season difficulty
                Gamma[i] = intersect(psi[i], np.vstack((psi[0][1:,:], np.zeros((1,I)))))
            else:
                Gamma[i] = intersect(psi[i], psi[i+1])
        verboseprint("Size", l, "done.", len(to_visit),
              "candidates to visit for size ", l+1, end='\r')
        l +=1
        to_visit = Lambda[:]

    verboseprint("End of MSGP_seasons")
    return res

def MSGP_seasons(data, minSup, closed=False, verbose=True):
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

    closed (boolean): if True return only closed gradual patterns with respect
    to seasons.

    Output:
    -------
    delta: a dictionary of all frequent seasonal gradual patterns
    '''
    # In case of verbose
    verboseprint = print if verbose else lambda *a, **k: None

    # Number of seasons, season length and number of attributes
    M, k, I = data.shape

    verboseprint("Beginning of MSGP_seasons:", M, k, I, minSup)

    Gamma = np.zeros((k, M, I), dtype=np.int8)
    verboseprint("Initialization...", end='\r')
    for stage in range(k):
        for cycle in range(M):
            if stage != k-1:
                Gamma[stage, cycle, :] = utils.get_gradual_items(data[cycle, [stage, stage+1], :])
            else:  # special treatment for last stage
                if cycle != M-1:  # use first stage of next season
                    Gamma[stage, cycle, :] = utils.get_gradual_items(np.vstack((data[cycle, stage, :], data[cycle+1, 0, :])))
                else:  # for last stage of last season, no data
                    Gamma[stage, cycle, :] = np.zeros(I)

    verboseprint("Starting BFS", end='\r')
    return MFCCP(Gamma, k, minSup, closed, verbose)

##################### Second algorithm #####################

def CFS(w, k, m):
    '''We suppose m != 1'''
    if not w:
        return set()

    sub = w[0]
    occ = dict([((start, 1), 0) for start in range(k)])
    occ[(sub, 1)] = 1

    l = 1
    for i in w[1:]:
        occ[(i, 1)] += 1
        if (i != (sub+l)%k):
            l = 1
            sub = i
        else:
            for j in range(l):
                if ((sub+j)%k, l-j+1) in occ:
                    occ[((sub+j)%k, l-j+1)] += 1
                else:
                    occ[((sub+j)%k, l-j+1)] = 1
            if l == k-1:
                sub = (sub+1)%k
            else:
                l += 1

    return set([season for season, occs in occ.items() if occs >= m])

def MFSP(Gamma, k, minSup):
    I = len(Gamma.keys())//2
    l = 1

    candidates = set([frozenset([i]) for i in range(2*I)])

    modulo = lambda d: d%k
    results = dict()
    while candidates and l <= I:
        psi = set()
        for p in candidates:
            seasons = CFS(list(map(modulo, Gamma[p])), k, minSup)

            if seasons:
                results[p] = seasons

                for p_prime in candidates:
                    g = p.union(p_prime)
                    if len(g) == l+1 and not g in Gamma:
                        Gamma[g] = sorted(set(Gamma[p]) & set(Gamma[p_prime]))
                        psi.add(g)

        candidates = psi
        l += 1

    return results

def MSGP_patterns(data, k, minSup, verbose=True):
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
    N, I = data.shape

    verboseprint = print if verbose else lambda *a, **k: None

    verboseprint("Beginning of MSGP_patterns:", N//k, k, I, minSup)
    Gamma = dict()
    for i in range(2*I):  # enumerate all gradual items
        Gamma[frozenset([i])] = list(utils.cover(i, data))

    return MFSP(Gamma, k, minSup)


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

    print("Beginning of Brute Force:", M, k, I, minSup)
    if format == 'season':
        results = dict()
        for start in range(k): # enumerate all seasons
            results[start] = dict()
            for length in range(1, k+1):
                # enumerate all gradual patterns
                grad_patterns = chain.from_iterable(combinations(range(2*I), l)
                                                    for l in range(1, 2*I+1))

                # Check each pattern over the current season
                for pattern in grad_patterns:
                    if utils.check_pattern((start, length), set(pattern),
                                           data.reshape((M,k,I)), minSup):
                        if results[start].get(length):
                            results[start][length].add(frozenset(pattern))
                        else:
                            results[start][length] = {frozenset(pattern)}

    else:
        results = dict()

        grad_patterns = chain.from_iterable(combinations(range(2*I), l)
                                                         for l in range(1, 2*I+1))
        for pattern in grad_patterns:  # enumerate all gradual patterns
            p = frozenset(pattern)
            for start in range(k):  # enumerate all seasons
                for length in range(1, k+1):
                    if utils.check_pattern((start, length), set(p),
                                           data.reshape((M, k, I)), minSup):
                        if results.get(p):
                            results[p].add((start, length))
                        else:
                            results[p] = {(start, length)}

    print("End of Brute Force")
    return results
