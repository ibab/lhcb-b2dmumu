import sys
import functools
import timeit
import logging

def prepare_sel(selections):
    return ' && '.join(map(lambda x: '(' + x + ')', selections))

def calc_tau(arr):
    from scipy.constants import c
    arr = arr['B_FD_OWNPV'] * arr['B_M'] / (arr['B_P'] * c * 10**3) * 10**12
    return arr

def binned_hist(ax, data, binedges, *args, **kwargs):
    #The dataset values are the bin centres
    x = (binedges[1:] + binedges[:-1]) / 2.0
    #The weights are the y-values of the input binned data
    weights = data
    return ax.hist(x, bins=binedges, weights=weights, *args, **kwargs)

def get_matching_variables(fname, tree, patterns):
    from fnmatch import fnmatch
    from root_numpy import list_branches

    branches = list_branches(fname, tree)

    selected = []

    for p in patterns:
        for b in branches:
            if fnmatch(b, p) and not b in selected:
                selected.append(b)
    return selected

def load_root(fname, tree=None, patterns=None, *kargs, **kwargs):
    """
    Loads a root file into a pandas DataFrame.
    Further *kargs and *kwargs are passed to root_numpy's root2array.

    >>> df = load_root('test.root', 'MyTree', patterns=['x_*', 'y_*'], selection='x_1 > 100')

    If the root file contains a branch called index, it will become the DataFrame's index.
    """
    from pandas import DataFrame
    from root_numpy import root2array, list_trees

    if tree == None:
        branches = list_trees(fname)
        if len(branches) == 1:
            tree = branches[0]
        else:
            raise ValueError('More than one tree found in {}'.format(fname))

    if not patterns:
        all_vars = None
    else:
        # index is always loaded if it exists
        patterns.append('index')
        all_vars = get_matching_variables(fname, tree, patterns)

    arr = root2array(fname, tree, all_vars, *kargs, **kwargs)
    if 'index' in arr.dtype.names:
        df = DataFrame.from_records(arr, index='index')
    else:
        df = DataFrame.from_records(arr)
    return df

def save_root(df, fname, tree_name, *kargs, **kwargs):
    from root_numpy import array2root
    logging.debug(df.columns)
    arr = df.to_records()
    array2root(arr, fname, tree_name, 'recreate', *kargs, **kwargs)

def time_job(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):

        result = [None]

        def sub_wrapper():
            result[0] = func(*args, **kwargs)

        timing = timeit.timeit(sub_wrapper, number=1)
        logging.warning('Execution of {} took {} secs'.format(func.__name__, timing))
        
        return result[0]
    return wrapper
