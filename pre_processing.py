'''
Utilitarian functions to pre-process data: missing values, entries, reformatting.
'''
import pandas as pd
import numpy as np
import csv

def stock_exchange_preproc(source, dest):
    '''Pre-process the stock exchange database.'''
    pass

def airQuality_preproc(source, dest):
    '''Process the AirQuality (https://archive.ics.uci.edu/ml/datasets/Air+Quality#)
    and save it as numpy binary file.
    source, dest are filepaths.'''
    # Data has a cycle length of 24, 389 seasons, 13 attributes
    data = pd.read_csv(source, delimiter=';', thousands=',')
    # Drop null value rows, and last incomplete season
    data.drop(range(9336,9471), inplace=True)
    # Drop date and time columns, as well as null value columns
    data.drop(columns=data.columns[[0,1,-2,-1]], inplace=True)
    np.save(dest, data.to_numpy().astype(np.float32))
