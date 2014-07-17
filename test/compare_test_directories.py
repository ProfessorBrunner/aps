#!/usr/bin/python
"""
Functions to judge the similarity of files using summary statistics. To
be used to compare the results of the shared memory and distributed
memory implementations.
"""
from os import listdir, system
from os.path import isfile, join, getmtime, isdir
from pandas import Series
from scipy.stats import kurtosis, tstd
from tabulate import tabulate
from math import sqrt
import numpy as np
import argparse
import matplotlib.pyplot as plt
#import matplotlib.colors.SymLogNorm as SymLogNorm
from SymLogNorm import SymLogNorm
import re

BINS = 605
TABLE_HEADERS = ['min', '25%', '50%', '75%', 'max']
IMPORTANT_FILES = [
    'signal.dat', 
    'iter_[0-9]*_C.dat',
    'iter_[0-9]*_fisher.dat',
    'iter_[0-9]*_window.dat']
IMPORTANT_REGEX_LIST = [re.compile(x) for x in IMPORTANT_FILES]

GRAPH_FILES = [
        'signal[0-9]{3}.dat',
        'iter_[0-9]*_fisher.dat', 
        'iter_[0-9]*_window.dat']
GRAPH_REGEX_LIST = [re.compile(x) for x in GRAPH_FILES]

def load_binary_file(file_name):
    """
    Load a binary file of doubles into a numpy array
    """
    return np.fromfile(file_name, dtype=np.dtype('f8'))

def matches_regex_list(str, regex_list):
    for regex in regex_list:
        if regex.match(str):
            return True
    return False

def direct_list_compare(expected, observed):
    print tabulate(zip(expected, observed),
            headers = ["Expected", "Observed"], floatfmt=".4f")
    
def print_numpy_with_format(matrix):
    for r in matrix:
        print ' '.join(['%3.10f' % x for x in r])

def raw_matrix(expected, observed):
    np.set_printoptions(threshold=np.nan, linewidth=8000)
    print "Expected"
    print_numpy_with_format(expected)
    print "Observed"
    print_numpy_with_format(observed)
    print "difference"
    print_numpy_with_format(expected - observed)

def get_file_error(expected_file, observed_file, reshape=False):
    """
    Load two binary data files and calculate and calculate errors.
    Args:
    expected_file: path to the "true" value binary file
    observed_file: path to the observed value binary file
    reshape: Boolean, if True attempt return results reshaped as square matrix
    """
    expected = load_binary_file(expected_file)
    observed = load_binary_file(observed_file)
    if not expected.size == observed.size:
        print "Size mismatch at {} and {}".format(expected_file, observed_file)
        return None
    if reshape:
        n = sqrt(len(expected))
        if n != int(n):
            print "Matrix is not square {}".format(expected_file, observed_file)
            return None
        expected = np.reshape(expected, (n, n))
        observed = np.reshape(observed, (n, n))

    error = expected - observed
    relative_error = error / expected

    return (expected, observed, error, relative_error)

def custom_compare_files(expected_path, observed_path, files):
    """
    This function is edited to do analysis of given files for debugging.
    """
    n = 0
    idx = range(n*BINS,(n+1)*BINS)
    for f in files:
        expected, observed, error, relative_error = \
            get_file_error(join(expected_path, f), join(observed_path, f), True)
        # raw_matrix(expected, observed)
        # # idx_exp = np.argsort(expected[0])
        # # idx_obs = np.argsort(observed[0])
        # # expected = expected[:, idx_exp]
        # # observed = observed[:, idx_obs]
        # # print tabulate(zip(expected[0], observed[0]),
        # #     headers = ["Expected", "Observed"])
        plt.pcolor(expected)
        plt.show()
        plt.pcolor(observed)
        plt.show()

        # for num in range(len(observed)):
        #     obs, exp = observed[num], expected[num]
        #     min_diff = min(abs(obs-exp), abs(obs+exp))
        #     if  min_diff > 1e-2:
        #         print "diff: {:8} exp: {:8} obs: {:8}".format(min_diff, expected[num], observed[num])
        #direct_list_compare(expected[0], observed[0])

def plot_heatmap(matrix, ax, label):
    vmax = np.max(matrix)
    vmin = np.min(matrix)
    vextreme = max(abs(vmin), vmax)
    k = kurtosis(matrix.flat)
    std = tstd(matrix.flat)
    #print "k: {} std: {}".format(k, std)
    args = {'vmax':vextreme,
            'vmin':-vextreme,
            'interpolation':'none',
            'aspect':'auto',
            'origin':'lower',
            'cmap':plt.get_cmap('RdBu')} #Spectral
    if k > 15:
        norm = SymLogNorm(std/3.0, vmin=-vextreme, vmax=vextreme)
        args['norm'] = norm
        label = "Symmetric log of " + label
    plt.imshow(matrix, **args)
    ax.set_title(label)
    ax.set_frame_on(False)
    plt.axis('off')
    # ax.grid(False)
    cb = plt.colorbar()
    if k > 15:
        ticks = np.linspace(0, 1, 9)
        tick_map = norm.inverse(ticks)
        cb.set_ticks(tick_map)
        cb.set_ticklabels(["{:.4g}".format(t) for t in tick_map])

def make_error_heatmap(expected_file, observed_file):
    try:
        expected, observed, error, relative_error = \
            get_file_error(expected_file, observed_file, reshape=True)
    except TypeError:
        return None
    fig = plt.figure()

    ax = plt.subplot(2, 1, 1)
    plot_heatmap(error, ax, "Error")

    ax = plt.subplot(2, 1, 2)
    plot_heatmap(relative_error, ax, "Relative Error")
    return fig

def make_error_boxplot(expected_files, observed_files, names):
    errors, relative_errors = [], []
    for expected_file, observed_file in zip(expected_files, observed_files):
        try:
            _, _, error, relative_error = \
                    get_file_error(expected_file, observed_file)
            errors.append(error)
            relative_errors.append(relative_error)
        except TypeError:
            return None

    fig = plt.figure()
    ax = plt.subplot(2, 1, 1)
    plt.boxplot(errors)
    plt.xticks(range(1, len(names)), names)

    ax = plt.subplot(2, 1, 2)
    plt.boxplot(relative_errors)
    plt.xticks(range(1, len(names)), names)

    return fig


def graph_directories(expected_path, observed_path, graph_type='box',
        file_regex_list=GRAPH_REGEX_LIST, output_dir=None):
    """
    Go through the files in each directory and graph comparison of pairs.
    """
    expected_files = [f for f in listdir(expected_path) if
            matches_regex_list(f, file_regex_list) and
            isfile(join(expected_path, f))]
    expected_files.sort(key=lambda x: x)
    observed_files = [f for f in listdir(observed_path) if
            matches_regex_list(f, file_regex_list) and
            isfile(join(observed_path, f))]

    available_files = [f for f in expected_files if f in observed_files]

    if graph_type == 'heat':
        for f in available_files:
            matrix_name = f.split('.')[0]
            print "Plotting {}".format(matrix_name)
            fig = make_error_heatmap(
                    join(expected_path, f), join(observed_path, f))
            if not fig:
                continue

            fig.suptitle(matrix_name)
            if output_dir:
                fig.savefig("{}/{}_heatmap".format(output_dir, matrix_name))
            else:
                fig.show()

    elif graph_type == 'box':
        groups, names = [], []
        for regex in file_regex_list:
            group = [f for f in available_files if regex.match(f)]
            if group:
                groups.append(group)
                names.append(group[0].split('.')[0])
        for group, matrix_name in zip(groups, names):
            print "Plotting {}".format(' '.join(group))
            fig = make_error_boxplot([join(expected_path, f) for f in group],
                    [join(observed_path, f) for f in group], group)
            if not fig:
                continue

            if output_dir:
                fig.savefig("{}/{}_boxplot".format(output_dir, matrix_name))
            else:
                fig.show()




#Also available: count, mean, std
def file_difference_summary_table(expected_file, observed_file):
    """
    Compare two files to check how similar they are.
    returns
    {expected:       [min, 25%, 50%, 75%, max],
     observed:       [min, 25%, 50%, 75%, max],
     error:          [min, 25%, 50%, 75%, max],
     relative_error: [min, 25%, 50%, 75%, max]}
    """
    try:
        expected, observed, error, relative_error = \
                get_file_error(expected_file, observed_file)
    except TypeError:
        return None

    # print "error>.1"
    # error3d = np.reshape(error, (BINS, BINS, 7), order='F')
    # idx = np.where(error3d>.05)
    # # row, col, bands = idx
    # for i, j, band in zip(*idx):
    #     print "i: {} j: {} band: {}".format(i,j,band)

    result = []
    for array_name in ('expected', 'observed', 'error', 'relative_error'):
        summary = Series(locals()[array_name]).describe()
        row = [array_name]
        row.extend([summary[header] for header in TABLE_HEADERS])
        result.append(row)

    return result

def compare_directories(expected_path, observed_path,
        file_regex_list=IMPORTANT_REGEX_LIST):
    """
    Go through the files in each directory and compare those that have
    matching names.

    Prints the resulting summary statistics.
    """
    expected_files = [f for f in listdir(expected_path) if
            matches_regex_list(f, file_regex_list) and
            isfile(join(expected_path, f))]
    expected_files.sort(key=lambda x: x)
    observed_files = [f for f in listdir(observed_path) if
            matches_regex_list(f, file_regex_list) and
            isfile(join(observed_path, f))]
    missing_files = []
    divider = '-'*80

    for f in expected_files:
        if f in observed_files:
            result = file_difference_summary_table(join(expected_path, f),
                    join(observed_path, f))
            print ""
            #print divider
            if result:
                print tabulate(result, headers=[f]+TABLE_HEADERS, 
                        tablefmt="orgtbl")
            observed_files.remove(f)
        else:
            missing_files.append(f)

    #print unmatched files
    if missing_files:
        print
        print "Following files in {} have no match in {}".format(expected_path,
                observed_path)
        print missing_files
    if observed_files:
        print
        print "Following files in {} have no match in {}".format(observed_path,
                expected_path)
        print observed_files


def main():
    """main method for using as a script"""
    #Command line parser
    parser = argparse.ArgumentParser(
        description='Compare binary files in two directories.')
    parser.add_argument("standard_dir", help="directory of accepted output.")
    parser.add_argument("testable_dir", help="directory of output to test.")
    parser.add_argument("-f", "--file", nargs='+', dest="specific_files",
        help="compare a specific file")
    parser.add_argument("-g", "--graph", dest="graph", action="store_true",
        help="display heatmap of important plots of data")
    parser.add_argument("-o", "--output", nargs=1, dest="output_dir",
        help="Output directory")
    args = parser.parse_args()
    if args.output_dir:
        args.output_dir = args.output_dir[0]
        if not isdir(args.output_dir):
            system('mkdir  -p ' + args.output_dir)
    print ""
    print ""

    if args.specific_files:
        custom_compare_files(args.standard_dir, args.testable_dir, 
            args.specific_files)
    elif args.graph:
        graph_directories(args.standard_dir, args.testable_dir, 
                output_dir=args.output_dir)
    else:
        compare_directories(args.standard_dir, args.testable_dir)


if __name__ == "__main__":
    main()
