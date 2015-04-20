import sys
import functools
import timeit
import logging
from pandas import DataFrame
from root_numpy import root2array, list_trees
from fnmatch import fnmatch
from root_numpy import list_branches

from luigi import LocalTarget

class RootTarget(LocalTarget):
    def __init__(self, path, tree):
        super(RootTarget, self).__init__(path)
        self.tree = tree

    def add_step(self, name):
        return RootTarget(self.path.replace('.root', '.{}.root'.format(name)), 'default')

def prepare_sel(selections):
    return ' & '.join(map(lambda x: '(' + x + ')', selections))

def calc_tau(arr):
    from scipy.constants import c
    arr = arr['B_FD_OWNPV'] * arr['B_M'] / (arr['B_P'] * c * 10**3) * 10**12
    return arr

def get_matching_variables(fname, tree, patterns):
    branches = list_branches(fname, tree)

    selected = []

    for p in patterns:
        for b in branches:
            if fnmatch(b, p) and not b in selected:
                selected.append(b)
    return selected

