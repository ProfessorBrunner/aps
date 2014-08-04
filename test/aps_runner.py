#!/usr/bin/python
"""
Create run scripts for performance and verification.
"""

from copy import deepcopy
from hashlib import sha224
import cPickle as pickle

EMAIL_ADDRESS="amwarren@email.arizona.edu"
#APS_DIR="/home/amwarren/aps"
APS_DIR="/home/alex/pop_uiuc/aps"
DATA_DIR="{}/data".format(APS_DIR)
OUTPUT_DIR="{}/output_distributed".format(DATA_DIR)
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
    cmd = "mpiexec -n {} ./distributed_memory/aps {}/{} {}/{} {}".format(
            run.mpi_nodes, DATA_DIR, run.fits_file, DATA_DIR, run.bands_file, run.name)
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
    'time':"00:30:00",
    'memory':"m24G",
    'cluster':"taub",
    'kl':False,
    'noise_model':'standard',
    'cross_correlation':False,
    'fits_file':"32_1000000_model_4.fits",
    'bands_file':"CL_32_model_4.bands",
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

runs = [aps_run(NUM_CORE_COMPARE)]
runs = cross_runs(runs, 'nodes', [1,2])
runs = cross_runs(runs, 'threads', [6,12])
runs = cross_runs(runs, 'mpi_nodes', [1,2])
runs = cross_runs(runs, 'repeat', [1,2,3])

submit_script = open("submit.bash", 'w')
submit_script.write('## Submission script\necho "## Abort script" > abort.bash\n')
for run in runs:
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

batch_name = create_batch_name(runs, "32_node_compare")
pickle_file = open("{}/{}.pkl".format(OUTPUT_DIR, batch_name), 'wb')
pickle.dump(runs, pickle_file)
