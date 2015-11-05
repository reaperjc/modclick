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

# import sys
# import subprocess
# import os
import time
# from utils import *
from functions import *
import Inputs
# import cPickle



######################
## PROGRAM COURSE ##
######################

if __name__ == '__main__':

    flags = sys.argv

    # Cheking inputs vars
    # inputs, superposition_pairs = prepareInputs(flags)
    inputs = Inputs.Inputs(flags)



    # saving tables info
    # d_table, d_table_inverse, d_pdbResnumTable = tables(inputs)

    # opening file for write outputs if user gave -o option
    if inputs.o:
        f_out = open(inputs.output, "w")
    else:
        f_out = open("results.output", "w")

    localtime = time.asctime(time.localtime(time.time()))
    f_out.write("Init time :\t%s\n" % localtime)




    # if list flag -l is given we add superposition pairs to superposition_pair list
    if '-l' in flags:
        add_superposition_pairs(inputs, superposition_pairs)

    results = []

    # first we generate all modfiles
    for mobile, target in superposition_pairs:
        if check(mobile["pdb"]) and check(target["pdb"]):
            # creating modified pdb files according the table
            modifyPDBAtoms(flags, mobile, d_table, d_pdbResnumTable)
            modifyPDBAtoms(flags, target, d_table, d_pdbResnumTable)

    if '-p' in flags:
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

        superposition_pair_lists = divide_lists(superposition_pairs, inputs.numCPUs)
        # ~ cPickle.dump(superposition_pair_lists,open("aers.data","w"))

        jobs = [(i, job_server.submit(click, (i, flags, d_table, d_pdbResnumTable, d_table_inverse,), \
                                      (leer, formatAtomName, check, modifyPDBAtoms, mergeClickResults,
                                       deleteModifiedFiles,), \
                                      ("subprocess", "time"))) for i in superposition_pair_lists]
        # ~ print [j() for i,j in jobs]
        for i, job in jobs:
            aux = job()
            if aux:
                results.extend(aux)

        job_server.print_stats()

    else:
        # running with just one core :(
        results = click(superposition_pairs, flags, d_table, d_pdbResnumTable, d_table_inverse)

    # Delete modified files
    for mobile, target in superposition_pairs:
        deleteModifiedFiles(flags, mobile, target, True)

    write_results(results, f_out)
    f_out.close()
    # if -i flag is used, now we restore the parameters.inp file
    if '-i' in flags:
        replaceParametersInp(inputs.parameters, restore=True)
