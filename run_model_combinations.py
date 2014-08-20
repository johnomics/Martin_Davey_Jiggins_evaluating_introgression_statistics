#!/usr/bin/env python3

# run_model_combinations.py
# Run a set of combined models using shared_ancestry_simulator.R based on a CSV file of model parameter combinations

# Written for "Evaluating statistics for the identification of introgressed loci"
# by Simon H. Martin, John W. Davey and Chris D. Jiggins
# Simon Martin: shm45@cam.ac.uk
# John Davey:   jd626@cam.ac.uk
# November-December 2013
# Revised February-May 2014, August 2014

import argparse
from math import floor
from datetime import datetime
from os import mkdir, chdir
from collections import defaultdict
from subprocess import call

parser = argparse.ArgumentParser(description='''Run a set of combined models using shared_ancestry_simulator.R.

    Takes a CSV file of model parameter combinations as input''')

parser.add_argument('-m', '--modelfile',  type=str, required=True)
parser.add_argument('-w', '--windows',  type=int, default=1000)
parser.add_argument('-t', '--threads', type=int, default=1)
parser.add_argument('-s', '--scalerate', type=float, default = 0.01)
parser.add_argument('-r', '--recombination', type=int, default = None)
parser.add_argument('-p', '--pathtosimulator', type=str, default="..")
parser.add_argument('-l', '--lengthofwindow', type=int, default = 5000)

args= parser.parse_args()

model_spec = """seqgen:
    m: HKY"""

model_spec += "\n    l: " + str(args.lengthofwindow)
model_spec += "\n    s: " + str(args.scalerate)
model_spec += """
ms:
    pops:
        - {name: 1, seqs: 8}
        - {name: 2, seqs: 8}
        - {name: 3, seqs: 8}
        - {name: 4, seqs: 8}
    opts:"""

if args.recombination != None:
    model_spec += "\n        recombination: " + str(args.recombination)

model_spec += """
        split:
            - {units: 3,   pops: [4,1]}
"""

def output_models(model_spec, model_type, params, log):
    pops = '0,0'
    if model_type is "Background":
        pops = '2,1'
    elif model_type is "Alternate":
        pops = '2,3'

    printpops = ''.join(pops.split(','))

    for model, model_windows in params.items():
        t123, t2x = model.split('/')

        model_t123 = "            - {units: " + t123 + ", pops: [3,1]}\n"
        model_txx =  "            - {units: " + t2x + ", pops: [" + pops + "]}\n"

        model_def = '{}_t123-{}_t{}-{}'.format(model_type, t123, printpops, t2x)
        try:
            model = open('{}.yml'.format(model_def), 'w')
            model.write(model_spec + model_t123 + model_txx)
            model.close()
        except:
            raise

        command = "{}/shared_ancestry_simulator.R -w {} -t {} -c {}.yml:1".format(
                    args.pathtosimulator, model_windows, args.threads, model_def)

        log.write('{}\t{}\n'.format(command,datetime.now()))
        try:
            call(command, shell=True)
        except:
            raise


def make_model_windows(args, log):
    try:
        models = open("../{}".format(args.modelfile))
    except:
        raise

    background = defaultdict(int)
    alternate = defaultdict(int)
    combinations = []

    header = models.readline().rstrip()
    
    for line in models:
        modeltype, bt123, bt21, at123, at23, altprop = line.rstrip().split(',')

        altprop = float(altprop)
        bg_windows = int(args.windows * (1-altprop))
        alt_windows = int(args.windows * altprop)
        combinations.append({'modeltype': modeltype,
                             'bt123': bt123, 'bt21': bt21,
                             'at123': at123, 'at23': at23,
                             'bg_windows': bg_windows,
                             'alt_windows':alt_windows})

        if bt123 != 'NA':
            pars = '{}/{}'.format(bt123,bt21)
            background[pars] = max(bg_windows, background[pars])
        if at123 != 'NA':
            pars = '{}/{}'.format(at123,at23)
            alternate[pars] = max(alt_windows, alternate[pars])

    output_models(model_spec, "Background", background, log)
    output_models(model_spec, "Alternate", alternate, log)

    return combinations, background, alternate


def write_windows(combf, model_filename, windows, startwin):
    modelf = open(model_filename)
    header = modelf.readline()
    if startwin == 0:
        combf.write(header)
    
    for win_i in range(0, windows):
        model = modelf.readline()
        block, blockNumber, blockName, *modelpars = model.split(',')
        blockNumber = int(blockNumber) + startwin
        
        outnum = '"{0}",{0},{0}'.format(blockNumber)
        
        combf.write(','.join([outnum] + modelpars))

    modelf.close()


def make_model_combinations(args, combinations, background, alternate):
    for comb in combinations:
        if comb['at123'] == 'NA': continue
            
        alt_pars = '{}/{}'.format(comb['at123'], comb['at23'])
        bg_pars = '{}/{}'.format(comb['bt123'], comb['bt21'])
        alt_spec = 'Alternate_t123-{}_t23-{}'.format(comb['at123'], comb['at23'])
        bg_spec = 'Background_t123-{}_t21-{}'.format(comb['bt123'], comb['bt21'])
            
        for seqs in 'ms', 'sg':
            combfilename = '{}.{}.{}.{}.csv'.format(alt_spec, bg_spec, args.windows, seqs)
            bg_filename = '{}.{}.{}.csv'.format(bg_spec, background[bg_pars], seqs)
            alt_filename = '{}.{}.{}.csv'.format(alt_spec, alternate[alt_pars], seqs)
            
            combf = open(combfilename,'w')
            write_windows(combf, bg_filename, comb['bg_windows'], 0)
            write_windows(combf, alt_filename, comb['alt_windows'], comb['bg_windows'])
            combf.close()

date = str(datetime.today()).split()[0]
outdir = "model_files_win{}_s{}_l{}_r{}_{}".format(args.windows, args.scalerate, args.lengthofwindow, args.recombination, date)

try:
    mkdir(outdir)
    chdir(outdir)
    log = open('{}.log'.format(outdir), 'w')
except:
    raise


combinations, background, alternate = make_model_windows(args, log)
make_model_combinations(args, combinations, background, alternate)