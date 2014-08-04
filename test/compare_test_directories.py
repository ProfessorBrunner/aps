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
from string import digits

SMALLEST_ERROR = 1e-14
BINS = 605
TABLE_HEADERS = ['min', '25%', '50%', '75%', 'max']
IMPORTANT_FILES = [
    # 'signal.dat',
    'C_iter_[0-9]*',
    'pre_window',
    'Y',
    'difference',
    'average',
    # # 'signal00[0-3]',
    'covariance_model_iter_',
    # 'kl_signal001',
    # # 'eigenvalues',
    # # 'eigenvectors',
    # # 'kl_noise',
    # # 'kl_overdensity',
    'fisher_iter_[0-9]*', 
    'window_iter_[0-9]*',
    ]
IMPORTANT_REGEX_LIST = [re.compile(x) for x in IMPORTANT_FILES]

GRAPH_FILES = [
    # 'signal[0-9]{3}',
    'inv_sqrt_fisher',
    'inv_fisher',
    'pre_window',
    'Y',
    'kl_signal[0-9]{3}',
    'average',
    # 'kl_noise',
    # 'kl_overdensity',
    'fisher_iter_[0-9]*', 
    'window_iter_[0-9]*',
    'C_iter_[0-9]*',
    'covariance_model_iter_[0-9]*',
    'difference',
    'preeigen',
    ]
# GRAPH_FILES=[
#     'kl_noise',
#     'kl_signal001',
#     'kl_overdensity',
# ]
GRAPH_REGEX_LIST = [re.compile(x) for x in GRAPH_FILES]

def format_name(name, strip_numbers=False):
    name = name.split('.')[0]
    name = name.replace('_', ' ')
    if strip_numbers:
        name = name.translate(None, digits)
        name = name.replace('iter', '')
        name = name.strip().capitalize()
    else:
        name = name.replace('iter', 'iteration')
    return name

def extract_number(string):
    return re.sub('\D', '', string)


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
        if n == int(n):
            n = int(n)
            expected = np.reshape(expected, (n, n), order='F')
            observed = np.reshape(observed, (n, n), order='F')
            #observed = expected.T

    # for row in xrange(0, n):
    #     for col in xrange(row+1, n):
    #         observed[col,row] = observed[row,col]

    error = expected - observed
    #error = np.where(abs(error) < SMALLEST_ERROR, error, 0)
    relative_error = error / expected

    return (expected, observed, error, relative_error)

THRESHOLD = .0001
def custom_compare_files(expected_path, observed_path, files):
    """
    This function is edited to do analysis of given files for debugging.
    """
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
        # comparison 
        # for col in xrange(len(expected[0])):
        #     diff = expected[:,col] - observed[:,col]
        #     combined = expected[:,col] + observed[:,col]
        for col_exp in xrange(len(expected[0])):
            distances = []
            neg_distances = []
            for col_obs in xrange(len(expected[0])):
                distances.append(np.linalg.norm(
                    expected[:,col_exp] - observed[:,col_obs]))
                neg_distances.append(np.linalg.norm(
                    expected[:,col_exp] + observed[:,col_obs]))
            #magnitude = np.linalg.norm(expected[:,col_exp])
            magnitude = .001 #this is the norm
            closest = min(distances)
            neg_closest = min(neg_distances)
            if closest < THRESHOLD and neg_closest < THRESHOLD:
                observed[:,col_exp] = 0
            elif closest < THRESHOLD:
                observed[:,col_exp] = .001
            elif neg_closest < THRESHOLD :
                observed[:,col_exp] = -.001
            #print_numpy_with_format([observed[:,col_exp],])
            #print "distance: {} {}".format(min(distances), magnitude)
            closest = min(closest, neg_closest)
            print "dist: %3.10f mag: %3.10f relative: %3.10f" \
                    % (closest, magnitude, closest/magnitude)

        ax = plt.subplot(1,1,1)
        fig = plot_heatmap(observed, ax, "observed")
        plt.show()
        # for num in xrange(len(observed)):
        #     obs, exp = observed[num], expected[num]
        #     min_diff = min(abs(obs-exp), abs(obs+exp))
        #     if  min_diff > 1e-2:
        #         print "diff: {:8} exp: {:8} obs: {:8}".format(min_diff, expected[num], observed[num])
        # direct_list_compare(expected[0], observed[0])

        # ax = plt.subplot(1,1,1)
        # fig = plot_heatmap(expected, ax, "expected")
        # plt.show()
        # ax = plt.subplot(1,1,1)
        # fig = plot_heatmap(observed, ax, "observed")
        # plt.show()


        # ax = plt.subplot(2,1,1)
        # fig = plot_heatmap(abs(error), ax, "abs error")
        # ax = plt.subplot(2,1,2)
        # fig = plot_heatmap(abs(relative_error), ax, "abs rel error")
        # plt.show()

def plot_heatmap(matrix, ax, label):
    vmax = np.max(matrix)
    vmin = np.min(matrix)
    vextreme = max(abs(vmin), vmax)
    k = kurtosis(matrix.flat)
    try:
        std = tstd(matrix.flat)
    except ZeroDivisionError:
        std = 0
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
    if len(matrix.shape) == 1:
        matrix = np.tile(matrix, (1, 2))
    plt.imshow(matrix, **args)
    ax.set_title(format_name(label))
    ax.set_frame_on(False)
    plt.axis('off')
    # ax.grid(False)
    ax.invert_yaxis()
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

    fig = plt.figure(figsize=(6,8))

    ax = plt.subplot(2, 1, 1)
    plot_heatmap(error, ax, "Error")

    ax = plt.subplot(2, 1, 2)
    plot_heatmap(relative_error, ax, "Relative Error")
    return fig

def make_compare_heatmap(expected_file, observed_file):
    try:
        expected, observed, error, relative_error = \
            get_file_error(expected_file, observed_file, reshape=True)
    except TypeError:
        return None

    fig = plt.figure(figsize=(12,10))
    ax = plt.subplot(2, 1, 1)
    plot_heatmap(expected, ax, "Expected")
    ax = plt.subplot(2, 1, 2)
    plot_heatmap(observed, ax, "Observed")
    return fig

def make_error_boxplot(expected_files, observed_files, names):
    #http://matplotlib.org/examples/pylab_examples/boxplot_demo2.html
    errors, relative_errors = [], []
    for expected_file, observed_file in zip(expected_files, observed_files):
        try:
            _, _, error, relative_error = \
                    get_file_error(expected_file, observed_file)
            errors.append(error)
            relative_errors.append(relative_error)
        except TypeError:
            return None

    fig = plt.figure(figsize=(6,4))
    ax = plt.subplot(2, 1, 1)
    plt.boxplot(errors)
    plt.xticks([])
    ax.set_title("Errors")
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
              alpha=0.5)

    ax = plt.subplot(2, 1, 2)
    plt.boxplot(relative_errors)
    ticks = [x + 1 for x in range(len(names))]
    names = [extract_number(name) for name in names]
    ax.set_title("Relative errors")
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
              alpha=0.5)

    return fig


def plot_band_powers(expected_path, observed_path, anafast_file, actual_file,
        output_dir=None):
    """Plot final bandpower result"""
    expected_files = [f for f in listdir(expected_path) if
            "C_iter_" in f and
            isfile(join(expected_path, f))]
    observed_files = [f for f in listdir(observed_path) if
            "C_iter_" in f and
            isfile(join(observed_path, f))]
    available_files = [f for f in expected_files if f in observed_files]

    _, band_center, _, _, anafast, _, _ = \
        np.loadtxt(anafast_file, dtype=float, unpack=True)
    anafast = np.sqrt(4*anafast)

    l, cl = np.loadtxt(actual_file, dtype=float, unpack=True)
    #cl = cl*2.*np.pi/(l*1.*(l+1.))
    cl = np.sqrt(cl)

    for current in available_files:
        name = format_name(current)
        expected = np.sqrt(4*load_binary_file(join(expected_path, current)))
        observed = np.sqrt(4*load_binary_file(join(observed_path, current)))

        #this is to compensate for the constant
        #added for convergence
        #anafast = anafast/1000

        fig = plt.figure()
        plt.plot(l, cl, label='True value')
        plt.plot(band_center, anafast, 'bo', label='Initial', markersize=5)
        plt.plot(band_center, expected, 'gx', label='Shared', markersize=10)
        plt.plot(band_center, observed, 'r+', label='Distributed', markersize=10)

        plt.yscale('log')
        plt.ylim(1e-3, 1e-0)
        plt.xlim(5, band_center.max())
        plt.xlabel(r'$\ell$', fontsize=18)
        plt.ylabel(r'$C_\ell$', fontsize=18)
        plt.legend(loc=0, numpoints=1)
        fig.suptitle("Bandpower Comparison {}".format(name), fontsize=24)
        if output_dir:
            fig.savefig("{}/C_anafast_comparison_{}".format(output_dir, name),
                        bbox_inches='tight', dpi=180)
            plt.close(fig)
        else:
            fig.show()
            fig.canvas.start_event_loop_default()


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

    if graph_type == 'heat' or graph_type == 'raw_heat':
        for f in available_files:
            matrix_name = f.split('.')[0]
            print "Plotting {}".format(matrix_name)
            if graph_type == 'heat':
                fig = make_error_heatmap(
                    join(expected_path, f), join(observed_path, f))
            elif graph_type == 'raw_heat':
                fig = make_compare_heatmap(
                    join(expected_path, f), join(observed_path, f))
            if not fig:
                continue

            fig.suptitle(format_name(matrix_name), fontsize=24, y=1.1)
            fig.canvas.set_window_title(matrix_name)
            if output_dir:
                fig.savefig("{}/{}_heatmap".format(output_dir, matrix_name),
                            bbox_inches='tight', dpi=180)
                plt.close(fig)
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
            if "iter" in matrix_name:
                plt.xlabel("Iteration")
            elif "signal" in matrix_name:
                plt.xlabel("Band")
            matrix_name = format_name(matrix_name, strip_numbers=True)
            fig.suptitle(matrix_name, fontsize=24, y=1.1)
            fig.canvas.set_window_title(matrix_name)
            if output_dir:
                fig.savefig("{}/{}_boxplot".format(output_dir, matrix_name), 
                        bbox_inches='tight', dpi=180)
                plt.close(fig)
            else:
                fig.show()
    if not output_dir:
        try:
            fig.canvas.start_event_loop_default()
        except:
            pass




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

    # observed = observed.flatten()
    # expected = expected.flatten()
    # error = error.flatten()
    # relative_error = expected.flatten()


    error = np.abs(error)
    relative_error = np.abs(relative_error)
    # print "error>.1"
    # error3d = np.reshape(error, (BINS, BINS, 7), order='F')
    # idx = np.where(error3d>.05)
    # # row, col, bands = idx
    # for i, j, band in zip(*idx):
    #     print "i: {} j: {} band: {}".format(i,j,band)

    result = []
    for array_name in ('expected', 'observed', 'error', 'relative_error'):
        #using locals is a hacky way of doing this
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
    parser.add_argument("--heat-plot", dest="plot_heat", action="store_true",
        help="display error heatmap of important plots of data")
    parser.add_argument("--heat-plot-compare", dest="plot_heat_compare",
        action="store_true", help="display heatmap of important plots of data")
    parser.add_argument("--box-plot", dest="plot_box", action="store_true",
        help="display error boxplot of important data")
    parser.add_argument("-c","-c-plot", nargs=2, dest="plot_c_file",
        help="Given anafast c values, compare with expected and observed")
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
    elif args.plot_heat:
        graph_directories(args.standard_dir, args.testable_dir, 
                output_dir=args.output_dir, graph_type='heat')
    elif args.plot_heat_compare:
        graph_directories(args.standard_dir, args.testable_dir, 
                output_dir=args.output_dir, graph_type='raw_heat')
    elif args.plot_box:
        graph_directories(args.standard_dir, args.testable_dir, 
                output_dir=args.output_dir, graph_type='box')
    elif args.plot_c_file:
        #args.plot_c_file = args.plot_c_file[0]
        plot_band_powers(args.standard_dir, args.testable_dir,
            args.plot_c_file[0], args.plot_c_file[1], output_dir=args.output_dir)
    else:
        compare_directories(args.standard_dir, args.testable_dir)

        


if __name__ == "__main__":
    main()
