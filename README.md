# Mining Seasonal Gradual Patterns

Author: Aymeric CÔME, November 2020 - April 2021

## About the project

These pieces of codes are the implementation of the algorithms I designed for extracting frequent seasonal gradual patterns from a temporal database.
The main reference is Extracting Seasonal Gradual Patterns from Temporal Sequence Data Using Periodic Patterns Mining, J Lonlac, A Doniec, M Lujak, SN Lecoeuche (https://arxiv.org/abs/2010.10289).
The algorithms have notably been inspired by the MPFPS BFS algorithm.
The pseudo-codes are given in the report I wrote; it hasn't been uploaded yet, but we're currently working on publishing it.

## Requirements

The code has been written and tested in Python 3.9.2, and uses Numpy 1.20.1.

## Documentation

In this section I explain a bit my code (functions, data structures).
We are dealing with cyclic temporal data sequences.
k is the cycle length, M is the number of seasons from the database, I is the number of attributes.
The timestamps and attributes are considered by their indices (fourth attribute corresponds to index 3).

### Data

The database used for experimentation part is the `AirQualityUCI.csv` one. The `pre_processing.py` script cleans it and save it in `AirQualityUCI.npy` as a 2D numpy array.
This dataset contains 389 cycles, each cycle contains 24 observations over 13 attributes.

### Seasons

This section is about the approach by seasons `MSGP_seasons`.

The input database should be a 3D numpy array: first dimension for cycles, second for stages inside cycles, last for attributes.
The graduality of an attribute is encoded by 1 if increasing, -1 if decreasing and 0 otherwise.
A transaction database collection (`Gamma`) is encoded as a 3D numpy array ((k, M, I)).
First dimension corresponds to a starting stage for the season considered (e.g. 2 when starting at the third stage of the cycle).
Second dimension corresponds to the cycles.
Third corresponds to the attributes.
For a season size l, `Gamma[i][j][m]` is the graduality of attribute `m` over the season of size l starting from stage `i`, in cycle `j`.
A gradual pattern is encoded as a (frozen)set: increasing attributes are specified by their indices and decreasing ones by their indices plus I (number of attributes. E.g. for I=6, {0, 2, 8} is the gradual pattern "attributes 0 and 2 are increasing and attribute 2 is decreasing".
The results are stored in a dictionary of dictionaries, mapping each starting stage and length (which defines a season) to the set of according frequent gradual patterns (res[start][k] = patterns).
There is a difficulty in the algorithm, when considering inter-season periods: indeed, to compute the corresponding gradual signature one has to fetch data from the following season. This must be considered at the initialization, and in the intersection.

To run the function over the Air Quality dataset, one simply has to load data `data = np.load('AirQualityUCI.npy)` and then `MSGP_seasons(data.reshape((389, 24, 13)), minSup)`.

### Gradual patterns

This section is about the approach by gradual patterns`'MSGP_patterns`.

The input database should be a 2D numpy array: first dimension for observations, second for attributes. The cycle length `k` is passed as an argument.
Gradual patterns are encoded as previously. A season is described as a tuple (starting_stage, length).
The time signatures are stored in a dictionary, as a list of indices.
'CFS' function returns the set of seasons from such a time signature.

Same procedure to use the algorithm, with `MSGP_patterns(data, k, minSup)`
