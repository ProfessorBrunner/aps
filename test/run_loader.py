#!/usr/bin/python
"""
Load run pickles and load output data for analysis
"""

from aps_runner import aps_run, APS_DIR, DATA_DIR, SOURCES_DIR, OUTPUT_DIR
import cPickle as pickle
import numpy as np
import matplotlib.pyplot as plt

INPUT_PICKLE="{}/32_node_compare_10d5b934602f6a17d736.pkl".format(OUTPUT_DIR)

def load_run(input_file):
    runs = pickle.load(open(input_file,'rb'))
    for run in runs:
        # print ""
        # print run
        # print "-"*30
        # print "nodes: %s threads: %s mpi_nodes: %s" % (run.nodes, run.threads, run.mpi_nodes)
        # print "-"*30
        run.load_result()
    return runs

def min_time(runs):
    print "times: "
    print [run.total_time for run in runs]
    return min(run.total_time for run in runs)

def max_mem(runs):
    print "mems: "
    print [run.torque_mem for run in runs]
    return max(run.torque_mem for run in runs)

def node_compare_matcher(nodes, threads, mpi_nodes):
    return lambda run : run.nodes == nodes and \
            run.threads == threads and run.mpi_nodes == mpi_nodes

def node_compare_lists(runs):
    times = []
    total_memory = []
    per_node_memory = []

    for num_nodes in [1, 2]:
        for num_cores in [1, 2]:
            for per_cores in [1, 2]:
                mpi_nodes = num_nodes*num_cores*per_cores
                print mpi_nodes
                threads = 6 * num_cores
                matcher = node_compare_matcher(num_nodes, threads, mpi_nodes)
                samples = filter(matcher, runs)
                times.append(min(run.total_time for run in samples))
                mem = max(run.torque_mem for run in samples) / float(2**20)
                total_memory.append(mem)
                per_node_memory.append(mem/float(num_nodes))
    print times
    print total_memory
    print per_node_memory

    N = len(times)
    ind = np.arange(N)
    width = .2

    fig = plt.figure(figsize=(8,6))
    ax = plt.subplot()
    fig.subplots_adjust(bottom=0.25, left=.15)
    memax = fig.add_axes(ax.get_position())
    ax.set_ylim(0, 1300)
    memax.set_ylim(0, 1300)
    memax.patch.set_visible(False)
    memax.yaxis.set_label_position('right')
    memax.yaxis.set_ticks_position('right')
    memax.spines['bottom'].set_position(('outward', 35))

    ax.yaxis.grid(True)
    ax.set_ylabel('Time (seconds)', fontsize=16)
    ax.get_xaxis().set_visible(False)
    memax.get_xaxis().set_visible(False)
    memax.set_ylabel('Memory (megabytes)', fontsize=16)
    memax.spines['bottom'].set_visible(False)


    rect_time = ax.bar(ind, times, width, color='r')
    rect_mem = memax.bar(ind+width, total_memory, width, color='b')
    rect_mem_per = memax.bar(ind+2*width, per_node_memory, width, color='c')
    #plt.axvline(x=3.8, ymin=-1, ymax = 1, linewidth=1, color='k')

    ax.legend( (rect_time[0], rect_mem[0], rect_mem_per[0]), 
            ('Runtime', 'Memory use', 'Memory per node'), 
            loc=9, prop={'size':16}, bbox_to_anchor=[.6, 1.])

    annotate_kwargs = {'xycoords':'axes fraction', 'textcoords':'offset points', 'va':'top', 'size':14}


    plt.annotate('Processes', (0,0), (-80, -15), **annotate_kwargs)
    interval = 440.0/N
    for i, val in zip(xrange(N), [1,2,2,4,2,4,4,8]):
        plt.annotate(str(val), (0,0), (15+i*interval, -15), **annotate_kwargs)

    plt.annotate('Cores', (0,0), (-80, -40), **annotate_kwargs)
    interval = 440.0/(N/2.0)
    for i, val in zip(xrange(4), [1,2,1,2]):
        plt.annotate(str(val), (0,0), (43+i*interval, -40), **annotate_kwargs)

    plt.annotate('Nodes', (0,0), (-80, -70), **annotate_kwargs)
    interval = 440.0/(N/4.0)
    for i, val in zip(xrange(2), [1,2]):
        plt.annotate(str(val), (0,0), (97+i*interval, -70), **annotate_kwargs)

    #plt.show()
    fig.savefig("run_compare.png", bbox_inches='tight', dpi=180)


def plot_cl(runs):
    l, cl = np.loadtxt("{}/{}".format(SOURCES_DIR, runs[0].source), dtype=float, unpack=True)
    cl = np.sqrt(cl)

    cls = [run.cl for run in runs]
    ls = [run.l for run in runs]

    fig = plt.figure(figsize=(8,6))

    plt.plot(l, cl, label='Model spectrum')
    plt.plot(ls[0], cls[0], 'r.', alpha=0.4, label="Estimated spectrum")
    for run_cl, run_l in zip(cls[1:], ls[1:]):
        plt.plot(run_l, run_cl, 'r.', alpha=0.4)
    plt.yscale('log')
    plt.ylim(1e-3, 1e-0)
    plt.xlim(5, run_l.max()+10)
    plt.xlabel(r'$\ell$', fontsize=18)
    plt.ylabel(r'$C_\ell$', fontsize=18)
    plt.legend(loc=0, numpoints=1)
    #plt.show()
    fig.savefig("cl.png", bbox_inches='tight', dpi=180)



if __name__ == "__main__":
    runs = load_run(INPUT_PICKLE)
    #plot_cl(runs)
    node_compare_lists(runs)

