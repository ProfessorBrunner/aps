#!/usr/bin/python
"""
Create run scripts for performance and verification.
"""

from generate_inputs_2 import aps_input
from copy import deepcopy
from hashlib import sha224
import cPickle as pickle
import numpy as np
import re

EMAIL_ADDRESS="amwarren@email.arizona.edu"
#APS_DIR="/home/amwarren/aps"
APS_DIR = "/home/alex/pop_uiuc/aps"
DATA_DIR = "{}/data".format(APS_DIR)
SOURCES_DIR = "{}/sources".format(DATA_DIR)
OUTPUT_DIR = "{}/output_distributed".format(DATA_DIR)
DIVIDER = "#" * 50

class aps_run:
    def __init__(self, kwargs):
        self.__dict__.update(kwargs)

    def update(self, d={}, **kwargs):
        self.__dict__.update(d)
        self.__dict__.update(kwargs)

    def duplicate(self):
        return deepcopy(self)

    def name_from_keys(self, keys, prefix=''):
        results = []
        if prefix:
            results.append(prefix)
        for key in keys:
            results.append(str(key) + "-" + str(self.__dict__[key]))
        results.append(self.create_unique_id())
        self.name = '_'.join(results)

    def create_unique_id(self):
        string = ""
        for key, val in self.__dict__.iteritems():
            if key != "hash_id" or key != "name":
                string += str(key) + str(val)
        self.hash_id = sha224(string).hexdigest()[0:10]
        return self.hash_id
    def __repr__(self):
        return self.name

    def load_result(self):
        self.time_results = {}
        self.mem_results = {}
        out_file = open("{}/out_{}".format(OUTPUT_DIR, self.name), 'r')
        re_time = re.compile(r'TIME:')
        re_malloc = re.compile(r'MALLOC:')
        re_resources = re.compile(r'Resources:')
        for line in out_file:
            line.strip()
            if re_time.match(line):
                self.load_result_time(line)
            elif re_malloc.match(line):
                self.load_result_malloc(line)
            elif re_resources.match(line):
                self.load_result_resources(line)

        self.malloc_max = 0
        self.malloc_max_id = ""

        for mem_id, mem_list in self.mem_results.iteritems():
            for mem in mem_list:
                if mem >  self.malloc_max:
                    self.malloc_max = mem
                    self.malloc_max_id = mem_id

        self.total_time = self.time_results['Total'][0]
        self.load_band_output()

        # print "%s: %d" % (self.malloc_max_id, self.malloc_max)
        # print "torque %d" % self.torque_mem
        # print "total time %f" % self.total_time

    def load_result_time(self, line):
        """load result for timing"""
        #print "T %s " % line
        line_split = line.split(' ')
        id_split = line_split[1].split('[')
        time = float(line_split[2])
        time_id = id_split[0]
        time_params = id_split[1]
        #print "%s    ~     %s " % (time_id, time)
        if time_id in self.time_results:
            self.time_results[time_id].append(time)
        else:
            self.time_results[time_id] = [time]

    def load_result_malloc(self, line):
        """load result for malloc"""
        #print "M %s " % line
        line_split = line.split(' ')
        id_split = line_split[1].split('.')
        mem_id = id_split[0]
        node = int(id_split[1])
        total_mem = int(line_split[3])

        if mem_id in self.mem_results:
            if node == 0:
                self.mem_results[mem_id].append(total_mem)
            else:
                self.mem_results[mem_id][-1] += total_mem
        else:
            self.mem_results[mem_id] = [total_mem]


    def load_result_resources(self, line):
        """load result for resource usage"""
        resources = line.split()
        resources = [x for x in resources if x]
        resources = resources[1].split(',')
        for resource in resources:
            resource_split = resource.split('=')
            key, val = resource_split[0], resource_split[1]
            if key == "mem":
                self.torque_mem = get_bytes(val)
            elif key == "vmem":
                self.torque_vmem = get_bytes(val)

    def load_band_output(self):
        file_name = "{}/C_{}.bands".format(OUTPUT_DIR, self.name)
        _, self.l, _, _, self.cl, _, _ = \
                np.loadtxt(file_name, dtype=float, unpack=True)
        neg = np.where(self.cl < 0)
        if neg[0].any():
            print "Negative at {}".format(self.name)
            print neg[0]
        self.cl = np.sqrt(4*self.cl)

def get_bytes(string):
    multiplier = 1
    if string[-1] == 'b':
        string = string[0:-1]
    if string[-1] == 'k':
        multiplier = 1024
        string = string[0:-1]
    return int(string)*multiplier

def cross_runs(runs, key, values):
    result = []
    for run in runs:
        result.append(run)
        run.__dict__[key] = values[0]
        for value in values[1:]:
            new_run = run.duplicate()
            new_run.__dict__[key] = value
            result.append(new_run)
    return result

def create_pbs(run):
    script = ["#!/bin/bash"]
    script.append("#")
    script.append("#PBS -N {}".format(run.name))
    script.append("#PBS -q secondary")
    script.append("#PBS -l naccesspolicy=singlejob")
    script.append("#PBS -l walltime={}".format(run.time))
    script.append("#PBS -l nodes={}:ppn={}:{}:{}".format(
            run.nodes, run.threads, run.memory, run.cluster))
    script.append("#PBS -m be")
    script.append("#PBS -M {}".format(EMAIL_ADDRESS))
    script.append("#PBS -j oe")
    script.append("#PBS -o {}/out_{}".format(OUTPUT_DIR, run.name))
    script.append("cd {}".format(APS_DIR))
    script.append("mkdir -p {}".format(OUTPUT_DIR))
    cmd = "mpiexec -verbose -n {} ./distributed_memory/aps {} {} {}".format(
            run.mpi_nodes, run.fits_file, run.bands_file, run.name)
    script.append('echo "{}"'.format(cmd))
    script.append(DIVIDER)
    script.append(cmd)
    return '\n'.join(script)

def create_batch_name(runs, prefix='batch'):
    names = [run.name for run in runs]
    unique = sha224(''.join(names)).hexdigest()[0:20]
    return "{}_{}".format(prefix, unique)



NUM_CORE_COMPARE = {
    'mpi_nodes':1,
    'nodes':1,
    'threads':12,
    'threads_per_core':6,
    'nside':32,
    'pixels':1024,
    'bands':40,
    'run':1,
    'time':"00:20:00",
    'memory':"m24G",
    'cluster':"taub",
    'kl':False,
    'noise_model':'standard',
    'cross_correlation':False,
    'fits_file':"32_1000000_model_4.fits",
    'bands_file':"CL_32_model_4.bands",
    'source':"CL_model0_0.30_0.40_zspec_799_fit2.dat",
    'name':"x",
}

STANDARD_ = {
    'nodes':1,
    'threads':12,
    'nside':64,
    'pixels':6836,
    'bands':10,
    'kl':False,
    'noise_model':'standard',
    'cross_correlation':False,
    'fits_file':"32_1000000_model_4.fits",
    'bands_file':"CL_32_model_4.bands",
    'name':"x",
}

def make_num_core_compare_batch():
    run = aps_run(NUM_CORE_COMPARE)
    aps_in = aps_input("{}/{}".format(SOURCES_DIR, run.source), run.nside)

    runs = [run]
    runs = cross_runs(runs, 'nodes', [1,2])
    runs = cross_runs(runs, 'threads', [6,12])
    runs = cross_runs(runs, 'mpi_nodes', [1,2])
    runs = cross_runs(runs, 'repeat', [1,2, 3])
    batch_name = create_batch_name(runs, "32_node_compare")


    submit_script = open("submit.bash", 'w')
    submit_script.write('## Submission script\necho "## Abort script" > abort.bash\n')
    for i, run in enumerate(runs):
        aps_in_temp = deepcopy(aps_in)
        aps_in_temp.patch_mask([int(x == i%12) for x in xrange(12)])
        run.fits_file = "{}/{}-{}.fits".format(DATA_DIR, batch_name, i)
        run.bands_file = "{}/{}-{}.bands".format(DATA_DIR, batch_name, i)

        aps_in_temp.write(run.fits_file, run.bands_file)

        run.pixels = aps_in_temp.pixels
        run.initial_cl = aps_in_temp.initial_cl

        run.cores = run.threads/run.threads_per_core * run.nodes
        run.mpi_nodes = run.mpi_nodes * run.cores
        run.name_from_keys(['nside', 'mpi_nodes', 'cores'], prefix="num_core_compare")

        pbs_file_name = "{}.pbs".format(run.name)
        pbs_file = open(pbs_file_name, 'w')

        qcmds = ['\n\n### {}'.format(run.name)]
        qcmds.append('jobname=`qsub {}.pbs`'.format(run.name))
        qcmds.append('echo "$jobname"')
        qcmds.append('echo "qdel $jobname" | cut -d\'.\' -f 1 >> abort.bash')
        submit_script.write('\n'.join(qcmds))

        pbs_file.write(create_pbs(run))

    pickle_file = open("{}/{}.pkl".format(OUTPUT_DIR, batch_name), 'wb')
    pickle.dump(runs, pickle_file)

if __name__ == "__main__":
    make_num_core_compare_batch()
