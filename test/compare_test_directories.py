#!/usr/bin/python
"""
"""
from os import listdir
from os.path import isfile, join, getsize
from numpy import savez, dtype , fromfile 
import argparse

def load_binary_file(file_name):
    """
    Load a binary file of doubles into a numpy array
    """
    pass

#http://docs.scipy.org/doc/numpy/reference/routines.io.html
#http://docs.scipy.org/doc/numpy/reference/generated/numpy.fromfile.html

def compare_files(expected_file, observed_file):
    """
    Compare two files to check how similar they are.
    """
    pass

def main():
    """main method for using as a script"""
    #Command line parser
    parser = argparse.ArgumentParser(
        description='Compare binary files in two directories.')
    parser.add_argument("standard_dir", help="directory of accepted output.")
    parser.add_argument("testable_dir", help="directory of output to test.")
    args = parser.parse_args()

def compare_directories(expected, observed):
    """
    """
    expected_files = [f for f in listdir(expected) if isfile(join(expected,f))]
    observed_files = [f for f in listdir(observed) if isfile(join(observed,f))]
    missing_files = []
    for f in expected_files:
        if f in observed_files:
            compare_files(join(expected_path, f), join(observed_path, f))
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


if __name__ == "__main__":
    main()
