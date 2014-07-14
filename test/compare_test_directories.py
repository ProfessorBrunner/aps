#!/usr/bin/python
"""
Functions to judge the similarity of files using summary statistics. To
be used to compare the results of the shared memory and distributed
memory implementations.
"""
from os import listdir
from os.path import isfile, join, getmtime
import numpy as np
from pandas import Series
import argparse
from tabulate import tabulate
import pylab as plt

BINS = 48

def load_binary_file(file_name):
    """
    Load a binary file of doubles into a numpy array
    """
    return np.fromfile(file_name, dtype=np.dtype('f8'))

#http://docs.scipy.org/doc/numpy/reference/routines.io.html
#http://docs.scipy.org/doc/numpy/reference/generated/numpy.fromfile.html

TABLE_HEADERS = ['min', '25%', '50%', '75%', 'max']
#Also available: count, mean, std
def compare_files(expected_file, observed_file):
    """
    Compare two files to check how similar they are.
    returns
    {expected:       [min, 25%, 50%, 75%, max],
     observed:       [min, 25%, 50%, 75%, max],
     error:          [min, 25%, 50%, 75%, max],
     relative_error: [min, 25%, 50%, 75%, max]}
    """
    expected = load_binary_file(expected_file)
    observed = load_binary_file(observed_file)
    if not expected.size == observed.size:
        print "Size mismatch at {} and {}".format(expected_file, observed_file)
        return None

    error = abs(expected-observed)

    relative_error = abs(error/expected)

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

def raw_table(expected, observed):
    print tabulate(zip(expected, observed),
            headers = ["Expected", "Observed"], floatfmt=".4f")

def raw_matrix(expected, observed):
    np.set_printoptions(threshold=np.nan, linewidth=8000)
    print "Expected"
    print expected
    print "Observed"
    print observed

def compare_files_table(expected_path, observed_path, files):
    """
    Compare two files and make a table summerizing thier differences
    """
    n =0
    idx = range(n*BINS,(n+1)*BINS)
    for f in files:
        expected = load_binary_file(join(expected_path, f)) #[idx]
        observed = load_binary_file(join(observed_path, f)) #[idx]
        # expected = np.reshape(expected, (BINS, BINS))
        # observed = np.reshape(observed, (BINS, BINS))
        expected = np.reshape(expected, (7, 7))
        observed = np.reshape(observed, (7, 7))
        raw_matrix(expected, observed)
        # # idx_exp = np.argsort(expected[0])
        # # idx_obs = np.argsort(observed[0])
        # # expected = expected[:, idx_exp]
        # # observed = observed[:, idx_obs]
        # # print tabulate(zip(expected[0], observed[0]),
        # #     headers = ["Expected", "Observed"])
        # plt.pcolor(expected)
        # plt.show()
        # plt.pcolor(observed)
        # plt.show()


        # for num in range(len(observed)):
        #     obs, exp = observed[num], expected[num]
        #     min_diff = min(abs(obs-exp), abs(obs+exp))
        #     if  min_diff > 1e-2:
        #         print "diff: {:8} exp: {:8} obs: {:8}".format(min_diff, expected[num], observed[num])
        raw_table(expected[0], observed[0])

def compare_directories(expected_path, observed_path):
    """
    Go through the files in each directory and compare those that have
    matching names.

    Prints the resulting summary statistics.
    """
    expected_files = [f for f in listdir(expected_path) if
            isfile(join(expected_path, f))]
    expected_files.sort(key=lambda x: getmtime(join(expected_path, x)))
    observed_files = [f for f in listdir(observed_path) if
            isfile(join(observed_path, f))]
    missing_files = []
    divider = '-'*80

    for f in expected_files:
        if f in observed_files:
            result = compare_files(join(expected_path, f),
                    join(observed_path, f))
            print ""
            #print divider
            if result:
                print tabulate(result, headers=[f]+TABLE_HEADERS, tablefmt="orgtbl")
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
    parser.add_argument("-f", "--file", nargs='+', dest="specific_file",
        help="compare a specific file")
    args = parser.parse_args()
    print ""
    print ""
    if args.specific_file:
        compare_files_table(args.standard_dir, args.testable_dir, 
            args.specific_file)
    else:
        compare_directories(args.standard_dir, args.testable_dir)


if __name__ == "__main__":
    main()
