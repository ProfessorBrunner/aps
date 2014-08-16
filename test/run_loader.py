#!/usr/bin/python
"""
Load run pickles and load output data for analysis
"""

from aps_runner import aps_run, APS_DIR, DATA_DIR, SOURCES_DIR, OUTPUT_DIR
import cPickle as pickle
import numpy as np
import matplotlib.pyplot as plt

INPUT_PICKLE="{}/32_node_compare_7a9f59c4db0592338363.pkl".format(OUTPUT_DIR)

def load_run(input_file):
    runs = pickle.load(open(input_file,'rb'))
    loaded_runs = []
    for run in runs:
        # print ""
        # print run
        # print "-"*30
        # print "nodes: %s threads: %s mpi_nodes: %s" % (run.nodes, run.threads, run.mpi_nodes)
        # print "-"*30
        try:
            run.load_result()
        except Exception as e:
            print("Error loading {}".format(run))
            print(e.message)
        else:
            loaded_runs.append(run)
    return loaded_runs

def min_time(runs):
    print "times: "
    print [run.total_time for run in runs]
    return min(run.total_time for run in runs)

def max_mem(runs):
    print "mems: "
    print [run.torque_mem for run in runs]
    return max(run.torque_mem for run in runs)

def node_compare_matcher(threads, mpi_nodes):
    return lambda run : run.threads == threads \
            and run.mpi_nodes == mpi_nodes

def node_compare_lists(runs):
    times = []
    total_memory = []
    per_node_memory = []

    for num_cores in [1, 2]:
        for mpi_nodes in [1,2,3,4,5,6,7,8]:
            threads = 6 * num_cores
            matcher = node_compare_matcher(threads, mpi_nodes)
            samples = filter(matcher, runs)
            print "Threads {} Processes: {} samples: {}".format(threads, mpi_nodes, len(samples))
            #print [run.total_time for run in samples]
            if len(samples) == 0:
                times.append(0.0)
                total_memory.append(0.0)
                continue
            times.append(min(run.total_time for run in samples))
            mem = max(run.torque_mem for run in samples) / float(2**20)
            total_memory.append(mem)
    print times
    print total_memory

    N = len(times)
    ind = np.arange(N)
    width = .2

    fig = plt.figure(figsize=(8,6))
    ax = plt.subplot()
    fig.subplots_adjust(bottom=0.25, left=.15)
    memax = fig.add_axes(ax.get_position())
    ax.set_ylim(0, 400)
    memax.set_ylim(0, 2400)
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
    #plt.axvline(x=3.8, ymin=-1, ymax = 1, linewidth=1, color='k')

    ax.legend( (rect_time[0], rect_mem[0]), 
            ('Runtime', 'Memory use'), 
            loc=9, prop={'size':16}, bbox_to_anchor=[.6, 1.])

    annotate_kwargs = {'xycoords':'axes fraction', 'textcoords':'offset points', 'va':'top', 'size':14}


    plt.annotate('Processes', (0,0), (-80, -15), **annotate_kwargs)
    interval = 440.0/N
    for i, val in zip(xrange(N), [1,2,2,4,2,4,4,8]):
        plt.annotate(str(val), (0,0), (15+i*interval, -15), **annotate_kwargs)

    plt.annotate('Cores', (0,0), (-80, -70), **annotate_kwargs)
    interval = 440.0/(N/4.0)
    for i, val in zip(xrange(2), [1,2]):
        plt.annotate(str(val), (0,0), (97+i*interval, -70), **annotate_kwargs)

    plt.show()
    #fig.savefig("run_compare.png", bbox_inches='tight', dpi=180)


def plot_cl(runs):
    l, cl = np.loadtxt("{}/{}".format(SOURCES_DIR, runs[0].source), dtype=float, unpack=True)
    cl = np.sqrt(cl)

    cls = [run.cl for run in runs]
    ls = [run.l for run in runs]

    fig = plt.figure(figsize=(8,6))

    plt.plot(l, cl, label='Model spectrum')
    plt.plot(ls[0], cls[0], 'g.', alpha=0.4, label="KL Estimated spectrum")
    for run_cl, run_l in zip(cls[1:], ls[1:]):
        plt.plot(run_l, run_cl, 'g.', alpha=0.4)
    plt.yscale('log')
    plt.ylim(1e-3, 1e-0)
    plt.xlim(5, run_l.max()+10)
    plt.xlabel(r'$\ell$', fontsize=18)
    plt.ylabel(r'$C_\ell$', fontsize=18)
    plt.legend(loc=0, numpoints=1)
    #plt.show()
    fig.savefig("kl_cl.png", bbox_inches='tight', dpi=180)



if __name__ == "__main__":
    runs = load_run(INPUT_PICKLE)
    #plot_cl(runs)
    node_compare_lists(runs)

