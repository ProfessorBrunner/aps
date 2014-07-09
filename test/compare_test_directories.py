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
    # error3d = np.reshape(error, (768, 768, 7), order='F')
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
        print "Following files in {} have no match in {}".format(expected_files,
                observed_files)
        for f in missing_files:
            print f
    if observed_files:
        print "Following files in {} have no match in {}".format(observed_files,
                expected_files)
        for f in observed_files:
            print f

def main():
    """main method for using as a script"""
    #Command line parser
    parser = argparse.ArgumentParser(
        description='Compare binary files in two directories.')
    parser.add_argument("standard_dir", help="directory of accepted output.")
    parser.add_argument("testable_dir", help="directory of output to test.")
    args = parser.parse_args()
    print ""
    print ""
    compare_directories(args.standard_dir, args.testable_dir)


if __name__ == "__main__":
    main()
