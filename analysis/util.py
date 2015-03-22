import sys
import functools
import timeit
import logging
from pandas import DataFrame
from root_numpy import root2array, list_trees
from fnmatch import fnmatch
from root_numpy import list_branches

def prepare_sel(selections):
    return ' & '.join(map(lambda x: '(' + x + ')', selections))

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
    branches = list_branches(fname, tree)

    selected = []

    for p in patterns:
        for b in branches:
            if fnmatch(b, p) and not b in selected:
                selected.append(b)
    return selected

def load_root(fname, tree=None, variables=None, ignore=None, *kargs, **kwargs):
    """
    Loads a root file into a pandas DataFrame.
    Further *kargs and *kwargs are passed to root_numpy's root2array.
    If the root file contains a branch called index, it will become the DataFrame's index.

    Params:
        fname - The filename of the root file
        tree - The name of the tree to load
        variables - A list of shell-patterns. Matching variables are loaded
        ignore - A list of shell-patterns. All matching variables are ignored (overriding the variables argument)

    Returns:
        A pandas DataFrame containing the loaded data

    Example:
        >>> df = load_root('test.root', 'MyTree', patterns=['x_*', 'y_*'], selection='x_1 > 100')

    """
    if tree == None:
        branches = list_trees(fname)
        if len(branches) == 1:
            tree = branches[0]
        else:
            raise ValueError('More than one tree found in {}'.format(fname))

    if not variables:
        all_vars = None
    else:
        # index is always loaded if it exists
        variables.append('index')
        all_vars = get_matching_variables(fname, tree, variables)

    if ignore:
        if not all_vars:
            all_vars = get_matching_variables(fname, tree, ['*'])

        ignored = get_matching_variables(fname, tree, ignore)
        if 'index' in ignored:
            raise ValueError('index variable is being ignored!')
        for var in ignored:
            all_vars.remove(var)

    arr = root2array(fname, tree, all_vars, *kargs, **kwargs)
    if 'index' in arr.dtype.names:
        df = DataFrame.from_records(arr, index='index')
    else:
        df = DataFrame.from_records(arr)
    return df

def save_root(df, fname, tree_name="default", *kargs, **kwargs):
    """
    Saves a pandas DataFrame as a root file.
    Further *kargs and *kwargs are passed to root_numpy's array2root.

    >>> df = DataFrame({'x': [1,2,3], 'y': [4,5,6]})
    >>> save_root('test.root', 'MyTree', df)
    
    The DataFrame index will be saved as an 'index' branch.
    """
    from root_numpy import array2root
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
