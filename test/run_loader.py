#!/usr/bin/python
"""
Load run pickles and load output data for analysis
"""

from aps_runner import aps_run, APS_DIR, DATA_DIR, SOURCES_DIR, OUTPUT_DIR
import cPickle as pickle
import numpy as np
import matplotlib.pyplot as plt

INPUT_PICKLE_2node="{}/32_node_compare_882796b23c43b9a3bafc.pkl".format(OUTPUT_DIR)
INPUT_PICKLE_1node="{}/32_node_compare_7a9f59c4db0592338363.pkl".format(OUTPUT_DIR)

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
    kl_times = []
    total_memory = []
    per_node_memory = []

    for num_cores in [1, 2]:
        for mpi_nodes in [1,2,3,4,5,6,7,8,9,10,11,12]:
            threads = 6 * num_cores
            matcher = node_compare_matcher(threads, mpi_nodes)
            samples = filter(matcher, runs)
            print "Threads {} Processes: {} samples: {}".format(threads, mpi_nodes, len(samples))
            #print [run.total_time for run in samples]
            if len(samples) == 0:
                times.append(0.0)
                kl_times.append(0.0)
                total_memory.append(0.0)
                continue
            kl_times.append(min(
                    run.time_results['KLCompression-eigensolve'][0] 
                    for run in samples
                    ))
            times.append(min(run.total_time for run in samples))
            mem = max(run.torque_mem for run in samples) / float(2**20)
            total_memory.append(mem)
            if not times[-1] == 0.0:
                print(kl_times[-1]/(times[-1]))
    print times
    print total_memory

    N = len(times)
    ind = np.arange(N)
    width = 1

    fig = plt.figure(figsize=(8,6))
    ax = plt.subplot()
    fig.subplots_adjust(bottom=0.25, left=.15)
    # memax = fig.add_axes(ax.get_position())
    ax.set_ylim(0, 1700)
    # memax.set_ylim(0, 2400)
    # memax.patch.set_visible(False)
    # memax.yaxis.set_label_position('right')
    # memax.yaxis.set_ticks_position('right')
    # memax.spines['bottom'].set_position(('outward', 35))

    ax.yaxis.grid(True)
    ax.set_ylabel('Time (seconds)', fontsize=16)
    ax.get_xaxis().set_visible(False)
    # memax.get_xaxis().set_visible(False)
    # memax.set_ylabel('Memory (megabytes)', fontsize=16)
    # memax.spines['bottom'].set_visible(False)


    rect_time = ax.bar(ind, times, width, color='r')
    rect_kl_time = ax.bar(ind, kl_times, width, color='g')
    # rect_mem = memax.bar(ind+width, total_memory, width, color='y')
    #plt.axvline(x=3.8, ymin=-1, ymax = 1, linewidth=1, color='k')

    ax.legend( (rect_time[0], rect_kl_time[0]), 
            ('Runtime', 'eigensolver time'), 
            loc=9, prop={'size':16}, bbox_to_anchor=[.3, 1.])

    annotate_kwargs = {'xycoords':'axes fraction', 'textcoords':'offset points', 'va':'top', 'size':14}

    plt.title("1 Nodes")
    plt.annotate('Processes', (0,0), (-80, -15), **annotate_kwargs)
    interval = 420.0/N
    for i, val in zip(xrange(N), [1,2,3,4,5,6,7,8,9,10,11,'',1,2,3,4,5,6,7,8,9,10,11,'']):
        plt.annotate(str(val), (0,0), (5+i*interval, -15), **annotate_kwargs)

    plt.annotate('Cores per node', (0,0), (-80, -70), **annotate_kwargs)
    interval = 420.0/(2)
    for i, val in zip(xrange(2), [1,2]):
        plt.annotate(str(val), (0,0), (97+i*interval, -70), **annotate_kwargs)

    plt.show()
    fig.savefig("1_node_comp_processes_times.png", bbox_inches='tight', dpi=180)


def aggregate_results(runs):
    times = []
    kl_times = []
    total_memory = []

    for num_cores in [1, 2]:
        for mpi_nodes in [1,2,3,4,5,6,7,8,9,10,11,12]:
            threads = 6 * num_cores
            matcher = node_compare_matcher(threads, mpi_nodes)
            samples = filter(matcher, runs)
            if len(samples) == 0:
                times.append(0.0)
                kl_times.append(0.0)
                total_memory.append(0.0)
                continue
            kl_times.append(min(
                    run.time_results['KLCompression-eigensolve'][0] 
                    for run in samples
                    ))
            times.append(min(run.total_time for run in samples))
            mem = max(run.torque_mem for run in samples) / float(2**20)
            total_memory.append(mem)

    return {'times':times,
            'kl_times':kl_times,
            'total_memory':total_memory}


def compare_node_memory(runs1, runs2):
    memory_1 = aggregate_results(runs1)['total_memory']
    memory_2 = aggregate_results(runs2)['total_memory']

    N = len(memory_1)
    ind = np.arange(N)
    width = .4

    fig = plt.figure(figsize=(8,6))
    ax = plt.subplot()
    fig.subplots_adjust(bottom=0.25, left=.15)
    ax.set_ylim(0, 1700)

    ax.yaxis.grid(True)
    ax.set_ylabel('Memory (megabytes)', fontsize=16)
    ax.get_xaxis().set_visible(False)

    rect_mem1 = ax.bar(ind, memory_1, width, color='r')
    rect_mem2 = ax.bar(ind+width, memory_2, width, color='g')

    ax.legend( (rect_mem1[0], rect_mem2[0]), 
            ('1 node', '2 nodes'), 
            loc=9, prop={'size':16}, bbox_to_anchor=[.3, 1.])

    annotate_kwargs = {'xycoords':'axes fraction', 'textcoords':'offset points', 'va':'top', 'size':14}

    plt.title("Memory comparison")
    plt.annotate('Processes', (0,0), (-80, -15), **annotate_kwargs)
    interval = 420.0/N
    for i, val in zip(xrange(N), [1,2,3,4,5,6,7,8,9,10,11,'',1,2,3,4,5,6,7,8,9,10,11,'']):
        plt.annotate(str(val), (0,0), (5+i*interval, -15), **annotate_kwargs)

    plt.annotate('Processors per node', (0,0), (-80, -70), **annotate_kwargs)
    interval = 420.0/(2)
    for i, val in zip(xrange(2), [1,2]):
        plt.annotate(str(val), (0,0), (97+i*interval, -70), **annotate_kwargs)

    plt.show()
    fig.savefig("node_comp_memory.png", bbox_inches='tight', dpi=180)

def compare_node_times(runs1, runs2):
    results_1 = aggregate_results(runs1)
    results_2 = aggregate_results(runs2)

    N = len(results_1['times'])
    ind = np.arange(N)
    width = .4

    fig = plt.figure(figsize=(8,6))
    ax = plt.subplot()
    fig.subplots_adjust(bottom=0.25, left=.15)
    ax.set_ylim(0, 1700)

    ax.yaxis.grid(True)
    ax.set_ylabel('Time (Seconds)', fontsize=16)
    ax.get_xaxis().set_visible(False)

    rect_time1 = ax.bar(ind, results_1['times'], width, color='#E8E28B')
    rect_kl_time1 = ax.bar(ind, results_1['kl_times'], width, color='#4C9645')
    
    rect_time2 = ax.bar(ind+width, results_2['times'], width, color='#C2BF9B')
    rect_kl_time2 = ax.bar(ind+width, results_2['kl_times'], width, color='#4F854A')

    ax.legend( (rect_time1[0], rect_kl_time1[0], rect_time2[0], rect_kl_time2[0]), 
            ('1 node time', 'eigensolver', '2 node time', 'eigensolver'), 
            loc=9, prop={'size':16}, bbox_to_anchor=[.2, 1.])

    annotate_kwargs = {'xycoords':'axes fraction', 'textcoords':'offset points', 'va':'top', 'size':14}

    plt.title("Time comparison")
    plt.annotate('Processes', (0,0), (-80, -15), **annotate_kwargs)
    interval = 420.0/N
    for i, val in zip(xrange(N), [1,2,3,4,5,6,7,8,9,10,11,'',1,2,3,4,5,6,7,8,9,10,11,'']):
        plt.annotate(str(val), (0,0), (5+i*interval, -15), **annotate_kwargs)

    plt.annotate('Processors per node', (0,0), (-80, -70), **annotate_kwargs)
    interval = 420.0/(2)
    for i, val in zip(xrange(2), [1,2]):
        plt.annotate(str(val), (0,0), (97+i*interval, -70), **annotate_kwargs)

    plt.show()
    fig.savefig("node_comp_time.png", bbox_inches='tight', dpi=180)

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
    #fig.savefig("kl_cl.png", bbox_inches='tight', dpi=180)



if __name__ == "__main__":
    runs1 = load_run(INPUT_PICKLE_1node)
    runs2 = load_run(INPUT_PICKLE_2node)
    compare_node_times(runs1, runs2)
    compare_node_memory(runs1, runs2)
    
    #plot_cl(runs)
    #node_compare_lists(runs)