"""
Date: 10-07-2012
Author: Jorge Cares
mail: jcaresgalvez@gmail.com
Description: The program is running as:
	python modclick.py pdb1.pdb pdb2.pdb -t tableFile -i ParametersFile
First, the programm creates modified versions of pdb files according the
table entered. Second, the programm run click with Parameters.inp modified
according the table file, and finnaly the programm merge the two click's
pdb files outputs separating the files in two MODELS
"""

import sys
# import subprocess
# import os
import time
# from utils import *
# from functions import *
import functions
import Inputs
import PDB

# import cPickle



######################
## PROGRAM COURSE ##
######################

if __name__ == '__main__':

    # Cheking inputs vars
    inputs = Inputs.Inputs(sys.argv)


    # opening file for write outputs if user gave -o option
    if inputs.o:
        f_out = open(inputs.output, "w")
    else:
        f_out = open("results.output", "w")

    localtime = time.asctime(time.localtime(time.time()))
    f_out.write("Init time :\t%s\n" % localtime)

    superposition_pairs = []  # list of PDB objects pairs that must be superposed

    # if -l flag is given, then we use superposition pairs given in the file
    if inputs.superpositions_list:
        functions.add_superposition_pairs(inputs, superposition_pairs)
    else:
        # if not, then use pdb1 and pdb2
        superposition_pairs.append((PDB.PDB(inputs.pdb1, inputs), PDB.PDB(inputs.pdb2, inputs)))


    results = []

    # first we generate all modfiles
    for mobile, target in superposition_pairs:
        if functions.check(mobile.pdb) and functions.check(target.pdb):
            # creating modified pdb files according the table
            # functions.modifyPDBAtoms(flags, mobile, inputs.d_table, inputs.d_pdbResnumTable)
            # functions.modifyPDBAtoms(flags, target, inputs.d_table, inputs.d_pdbResnumTable)
            mobile.modifyPDBAtoms(inputs)
            target.modifyPDBAtoms(inputs)
            # functions.modifyPDBAtoms2(inputs, mobile, fullpdb=None, pdbname=None,)
            # functions.modifyPDBAtoms2(inputs, target, fullpdb=None, pdbname=None,)

    # if -p flag is given then try parallel python
    if inputs.p:
        # with this flag we try to run parallel python
        try:
            import pp
        except:
            print "No pp module in python library\n"
            sys.exit(1)

        ppservers = ()

        if inputs.numCPUs:
            job_server = pp.Server(ncpus=inputs.numCPUs, ppservers=ppservers)
        else:
            job_server = pp.Server(ppservers=ppservers)
            inputs.numCPUs = job_server.get_ncpus()

        # now divide superpositions pair list to try to pass same jobs quantities to each core
        divided_superposition_pair_lists = functions.divide_lists(superposition_pairs, inputs.numCPUs)
        # ~ cPickle.dump(superposition_pair_lists,open("aers.data","w"))

        # list comprehension of parallel python jobs
        jobs = [(
                    i_superposition_pairs,
                    job_server.submit(
                        functions.click,
                        (i_superposition_pairs, inputs),
                        (
                            functions.leer,
                            functions.formatAtomName,
                            functions.check,
                            functions.modifyPDBAtoms,
                            functions.mergeClickResults,
                            functions.deleteModifiedFiles,
                        ),
                        ("subprocess", "time")))

                for i_superposition_pairs in divided_superposition_pair_lists
                ]

        # Run jobs
        # ~ print [j() for i,j in jobs]
        for i, job in jobs:
            aux = job()
            if aux:
                results.extend(aux)

        job_server.print_stats()

    else:
        # running with just one core :(
        results = functions.click(superposition_pairs, inputs)
        print results
        print "caca"

    # Delete modified files
    # for mobile, target in superposition_pairs:
    #     functions.deleteModifiedFiles(inputs, mobile, target, True)

    # wite results
    functions.write_results(results, f_out)
    f_out.close()


    # if -i flag is used, restore the parameters.inp file
    if inputs.i:
        inputs.replaceParametersInp(inputs.parameters, restore=True)
