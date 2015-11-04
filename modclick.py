'''
Date: 10-07-2012
Author: Jorge Cares
mail: jcaresgalvez@gmail.com
Description: The program is running as:
	python modclick.py pdb1.pdb pdb2.pdb -t tableFile -i ParametersFile
First, the programm creates modified versions of pdb files according the
table entered. Second, the programm run click with Parameters.inp modified
according the table file, and finnaly the programm merge the two click's 
pdb files outputs separating the files in two MODELS'''

import sys
import subprocess
import os
import time
import cPickle


#############
## Classes ##
#############

class Inputs:
    '''Inputs class creates a object with command line entry inputs'''

    def __init__(self):
        self.parameters = None
        self.table = None
        self.pdb_resnum_table = None
        self.superpositions_list = None
        self.numCPUs = None
        self.output = None


###############
## FUNCTIONS ##
###############

def PDB(pdb):
    '''Antes una clase, ahora una funcion que crea un diccionario'''
    d = {}
    d["pdb"] = os.path.abspath(pdb)
    abspathPDB = os.path.abspath(pdb)
    partes = abspathPDB.split("/")
    d["name"] = partes[-1][:-4]
    d["path"] = "/".join(partes[:-1] + ["", ])
    return d


def tables(inputs):
    '''Open the "table" file with new atomtypes specifications and save
	the info into a dicctionary "d_table". A inverse dictionary "d_table_inverse"
	is also created to restore PDB files to IUPAC nomenclature.'''

    d_table = {}  # dictionary with {atomName:newAtomType}
    d_table_inverse = {}  # dictionary inverse to d_table
    d_pdbResnumTable = {}  # dictionary with {PDB:{CHAIN:{resnum:None,},},}
    for i, line in enumerate(open(inputs.table)):
        if i > 0 and len(line.strip()) > 0:
            aux = line.split(":")
            if len(aux) == 3:
                if not d_table.has_key(aux[0].strip()):
                    d_table[aux[0].strip()] = {aux[1].strip(): aux[2].strip(), }
                    d_table_inverse[aux[0].strip()] = {aux[2].strip(): aux[1].strip(), }
                else:
                    d_table[aux[0].strip()][aux[1].strip()] = aux[2].strip()
                    d_table_inverse[aux[0].strip()][aux[2].strip()] = aux[1].strip()

    # now if -f flag is choosen
    if inputs.pdb_resnum_table:
        f = open(inputs.pdb_resnum_table)
        l = f.readlines()
        f.close()
        for i in l[1:]:
            aux = i.split()
            pdb = aux[0].strip()
            chain = aux[1].strip()
            number = aux[2].strip()
            if not d_pdbResnumTable.has_key(pdb):
                d_pdbResnumTable[aux[0]] = {}
            if not d_pdbResnumTable[aux[0]].has_key(chain):
                d_pdbResnumTable[pdb][chain] = {}
            d_pdbResnumTable[pdb][chain][number] = None
    else:
        d_pdbResnumTable = None

    return d_table, d_table_inverse, d_pdbResnumTable


def formatAtomName(atomName):
    '''Return the atomname field with the correct PDB format'''
    if len(atomName) == 1:
        return " %s  " % atomName
    elif len(atomName) == 2:
        return " %s " % atomName
    elif len(atomName) == 3:
        return " %s" % atomName
    elif len(atomName) == 4:
        return "%s" % atomName


def modifyPDBAtoms(flags, pdb_obj, d_table, d_pdbResnumTable, fullpdb=None, pdbname=None, clickOutput=False):
    '''Take a pdb file as input and returns a new PDB file with atomtypes
	modifications. The new file has the "Mod" word before .pdb extension.
	If clickOutput flag is True, the pdb file will be modified with normal 
	IUPAC names for a standard PDB file'''
    d_keywordRecord = {"ATOM": None, "HETATM": None, "ANISOU": None, "SIGUIJ": None}

    if fullpdb:
        pdb = fullpdb
    else:
        pdb = pdb_obj["pdb"]
    if pdbname:
        name_pdb = pdbname
    else:
        name_pdb = pdb_obj["name"]

    l_newPDB = []
    for line in open(pdb):
        if line[:6].strip() in d_keywordRecord:
            resname = line[17:20].strip()
            atomname = line[12:16].strip()
            chain = line[21]
            resnum = line[22:26].strip()

            # first we can modify a line just if the table has the right resname and atomname
            if d_table.has_key(resname):
                if d_table[resname].has_key(atomname):
                    if '-f' in flags:
                        # if -f flag is given, we use only the line with residue specified in d_pdbResnumTable
                        if name_pdb in d_pdbResnumTable and \
                                        chain in d_pdbResnumTable[name_pdb] and \
                                        resnum in d_pdbResnumTable[name_pdb][chain]:
                            l_newPDB.append(line[:12] + formatAtomName(d_table[resname][atomname]) + line[16:])
                        else:
                            # ~ l_newPDB.append(line[:12] + formatAtomName(d_table[resname][atomname]) + line[16:])
                            l_newPDB.append(line)
                    else:
                        l_newPDB.append(line[:12] + formatAtomName(d_table[resname][atomname]) + line[16:])
                else:
                    l_newPDB.append(line)
            else:
                l_newPDB.append(line)
        else:
            l_newPDB.append(line)
    if not clickOutput:
        f = open("%sMod.pdb" % pdb[:-4], "w")
    else:
        f = open("%s.pdb" % pdb[:-4], "w")
    for i in l_newPDB:
        f.write(i)
    f.close()


def modifyClickParameters():
    '''modify the Parameters.inp file according the table information'''
    # Creating new typeAtom line
    newTypeAtomLine = "typeAtom=\""
    typeatoms = []
    for i in d_table:
        typeatoms.extend(d_table[i].values())
    daux = {i: None for i in typeatoms}
    typeatoms = daux.keys()
    typeatoms.sort()
    for i, j in enumerate(typeatoms):
        if i == len(typeatoms) - 1:
            newTypeAtomLine += formatAtomName(j)
        else:
            newTypeAtomLine += "%s ," % formatAtomName(j)
    newTypeAtomLine += " \"\n"

    # Creating new line list for the new Parameters.inp file
    l_newLines = []
    for i, line in enumerate(open("Parameters.inp")):
        if i == 0:
            l_newLines.append(newTypeAtomLine)
        else:
            l_newLines.append(line)

    # Writing new Parameters.inp file
    f = open("Parameters.inp", "w")
    for i in l_newLines:
        f.write(i)
    f.close()


def deleteModifiedFiles(flags, mobile, target, deleteModfiles=False):
    '''This function deletes temp modified pdb files and the two traslated
	pdbs generated by click'''

    # deleting Mod pdb files created
    # ~ target["pdb"][:-4])],stdout=subprocess.PIPE ,shell = True, close_fds=True)
    if deleteModfiles:
        # ~ proc = subprocess.call(["rm -f %sMod.pdb %sMod.pdb "%(
        # ~ mobile["pdb"][:-4],
        # ~ target["pdb"][:-4])],
        # ~ shell = True)
        try:
            os.unlink("%sMod.pdb" % (mobile["pdb"][:-4]))
            os.unlink("%sMod.pdb" % (target["pdb"][:-4]))
        except:
            pass


    else:

        # ~ if not '-s' in flags:
        # ~ proc = subprocess.call(["rm -f %sMod-%sMod.1.pdb %sMod-%sMod.1.pdb "%(mobile["pdb"][:-4],
        # ~ target["name"],
        # ~ target["pdb"][:-4],
        # ~ mobile["name"])],
        # ~ shell = True)
        try:
            os.unlink("%sMod-%sMod.1.pdb" % (mobile["pdb"][:-4], target["name"]))
            os.unlink("%sMod-%sMod.1.pdb" % (target["pdb"][:-4], mobile["name"]))
        except:
            pass

        if not '-m' in flags:
            # ~ proc = subprocess.call(["rm -f %s%sMod-%sMod.pdb.1.clique"%(mobile["path"],
            # ~ mobile["name"],
            # ~ target["name"])],
            # ~ shell = True)
            try:
                os.unlink("%s%sMod-%sMod.pdb.1.clique" % (mobile["path"],
                                                          mobile["name"],
                                                          target["name"]))
                os.unlink("%s%s_%s_Merged.pdb" % (mobile["path"], mobile["name"], target["name"]))
            except:
                pass


def mergeClickResults(flags, mobile, target):
    '''Merge click results in a single file. The files are separated with
	MODEL record'''

    try:
        f2 = open("%s%sMod-%sMod.1.pdb" % (mobile["path"], mobile["name"], target["name"]))

        # ~ f3 = open("%s%sMod-%sMod.1.pdb"%(target["path"], target["name"], mobile["name"]))
        f3 = open("%s%sMod.pdb" % (target["path"], target["name"]))

    except:
        deleteModifiedFiles(flags, mobile, target)
        # ~ sys.exit(1)
        return 0

    f = open("%s%s_%s_Merged.pdb" % (mobile["path"], mobile["name"], target["name"]), "w")
    f.write("%s\n" % ("MODEL        1".ljust(80)))

    for i in f2.readlines():
        f.write(i)
    f.write("%s\n" % "ENDMDL".ljust(80))
    f.write("%s\n" % ("MODEL        2".ljust(80)))

    for i in f3.readlines():
        f.write(i)
    f.write("ENDMDL".ljust(80))

    f.close()
    f2.close()
    f3.close()


def check(input_file):
    '''Checks if inputs files exist'''
    # we check only if a file name was given
    if input_file:
        try:
            f = open(input_file)
            f.close()
            return True
        except:
            return False
    return True


def prepareInputs(inputs):
    '''Prepare all files according to flags inputs'''
    nom_atomos = {'A': ['N1', 'C2', 'N3', 'C4', 'C5', 'C6', 'N7', 'C8', 'N9', 'N6'], \
                  'G': ['N1', 'C2', 'N3', 'C4', 'C5', 'C6', 'N7', 'C8', 'N9', 'O6', 'N2'], \
                  'C': ['N1', 'C2', 'N3', 'C4', 'C5', 'C6', 'O2', 'N4'], \
                  'T': ['N1', 'C2', 'N3', 'C4', 'C5', 'C6', 'O2', 'O4', 'C7'],
				  # notar que el C7 es el metilo en el C5 \
                  'U': ['N1', 'C2', 'N3', 'C4', 'C5', 'C6', 'O2', 'O4']}
    # creating inputs object
    inputs = Inputs()
    superposition_pairs = []  # list of all superposition pairs
    par = []  # var use for save pdb inputs

    for j, i in enumerate(sys.argv):
        # check the pdb files
        if i.strip()[-4:] == '.pdb':
            par.append(i.strip())

        # checking parameters file
        if i.strip() == '-i':
            try:
                inputs.parameters = sys.argv[j + 1]
            except:
                print "No Parameters file indicated"
                errorSentence()

        # checking table file
        if i.strip() == '-t':
            try:
                inputs.table = sys.argv[j + 1]
            except:
                print "No table file indicated"
                errorSentence()

        # checking pdb_resnum_table
        if i.strip() == '-f':
            try:
                inputs.pdb_resnum_table = sys.argv[j + 1]
            except:
                print "No PDB_resnum_table file indicated"
                errorSentence()

        # checking superposition list
        if i.strip() == '-l':
            try:
                inputs.superpositions_list = sys.argv[j + 1]
            except:
                print "No superpositions list file indicated"
                errorSentence()

        # checking parallel python nCPUs
        if i.strip() == '-p':
            try:
                if sys.argv[j + 1].isdigit():
                    inputs.numCPUs = int(sys.argv[j + 1])
            except:
                pass

        # checking if output file name is given
        if i.strip() == '-o':
            try:
                inputs.output = sys.argv[j + 1].strip()
            except:
                pass

    # cheking if the files exists

    if not len(par) in (0, 2):
        print "two pdf files are necessary as minimum or a list of superposition pairs."
        errorSentence()

    if par and par[0] and par[1]:
        # if user gives two pdb files
        if not check(par[0]):
            print "%s does'n exist, please try a real pdb file." % par[0]
            errorSentence()

        if not check(par[1]):
            print "%s does'n exist, please try a real pdb file." % par[1]
            errorSentence()
        # if this pdbs exist, we create PDB objects
        superposition_pairs.append((PDB(par[0]), PDB(par[1])))
    else:
        # if no list of superposition is given
        if not '-l' in sys.argv:
            print "Click needs two pdb files."
            errorSentence()

    # if check
    if not check(inputs.parameters):
        print "%s does'n exist, please try a real parameters file." % inputs.parameters
        errorSentence()

    if not inputs.table:
        print "modclick needs a table file."
        errorSentence()
    else:
        if not check(inputs.table):
            print "%s does'n exist, please try a real table file." % inputs.table
            errorSentence()

    if not check(inputs.pdb_resnum_table):
        print "%s does'n exist, please try a real PDB_resnum_table file." % inputs.pdb_resnum_table
        errorSentence()

    if not check(inputs.superpositions_list):
        print "%s does'n exist, please try a real superposition list." % inputs.superpositions_list
        errorSentence()

    # checking Parameters.inp file
    if not check("Parameters.inp"):
        print "Parameters.inp does'n exist."
        sys.exit(1)

    # if everething goes fine
    return inputs, superposition_pairs


def errorSentence():
    # ~ print "\nTry one of this two options:\n1)python modclick pdb1.pdb pdb2.pdb -t tableFile"
    # ~ print "2)python modclick pdb1.pdb pdb2.pdb -t tableFile -i ParametersFile\n"
    helpInfo()
    sys.exit(1)


def replaceParametersInp(newParametersFile, restore=False):
    '''Realize a backup of Parameters.inp file, replace it with a new
	parameters file. If restore flag is True Parameters.inp backup is
	restored'''

    # creating a Parameters.inp.backup file
    if not restore:
        f = open("Parameters.inp.backup", "w")
        for i in open("Parameters.inp"):
            f.write(i)
        f.close()

        # replace Parameters.inp with the new file
        f = open("Parameters.inp", "w")
        for i in open(newParametersFile):
            f.write(i)
        f.close()
    else:
        f = open("Parameters.inp", "w")
        for i in open("Parameters.inp.backup"):
            f.write(i)
        f.close()
        proc = subprocess.Popen(["rm -f Parameters.inp.backup"], stdout=subprocess.PIPE, shell=True)
        l2 = proc.communicate()


def helpInfo():
    print """
DEPENDENCIES

	This python script needs the next files in the same script folder:
		*Click binay
		*Parameters.inp
		*table file specifing new atomtypes
		*optional: Another Parameters file
		*optional: table file with PDB, chain, resnum specification
		*optional: list with pdb1,pdb2 superposition pairs
	Additionally you need parallel python module to run several click
	instances using	more than one core in a multicore system.

MINIMAL USAGE: 

	python modclick.py mobile.pdb target.pdb -t tableFile

OPTIONS:

	-t [file] modclick realize superpositions with the new atomnames specified
		in the table [file] with resname:atomname:new_atomname specifitacions.
		modclick changes the typeAtom record of Parameters.inp file by the atom
		types of this table automatically

	-i [file] modclick uses the [file] parameters to run click superpositions

	-f [file] modclick uses the [file] table with PDB:chain:resnumber fields that
		shows only the residues wich you want to change atom types

	-l [file] modclick runs click over a pair list of PDB files separated by
		comma ',' from [file]. For all those PDB files the same table of atom
		types will be applied and click will run with the same parameters
	
	-m	modclick mantains *.clique file. Merged pdb file its mantained too except if
		-s flag is chosen
		
	-k  kept transformed pdb files created by click except if -s flag is chosen
	
	-s	modclick runs click with -s 0 flag. With this flag click doesn't create
		transformation pdb files therefore Merged pdb file is not created
	
	-p [nCPUs] modclick try to run superpositions with parallel python module using
		a given number of cores. If no number is specified, or 0 number is given,
		modclick uses the maximun machine available cores
	
	-o [file] modclick writes summary [file] with pdb1,pdb2,RMSD records. If not
		[file] name is given or this flag is not used, results.out file will be 
		written

EXAMPLE:
	python modclick.py -t MyAtomtypesTable -l MyList -s -p 4 -o MyOutput
	this line runs click on each superposition pair in MyList with MyParameters 
	Parameters using MyAtomtypesTable using 4 cores. Restults (RMSD) will be written in MyOutput file

OUTPUTS:

	mobile_target_Merged.pdb:	this file contain the 2 click's outputs merged
		in a single pdb file separating with MODEL records. MODEL 1 contain
		coordenates for mobile superposed on target, and MODEL 2 the target
		coordenates

	mobileMod-targetMod.pdb.1.clique:	It's the same click's output but the
		"Mod" word is added because click was running with modified pdb files.
		If both pdbs are in different folders, this file will be written in the
		folder of first pdb file
	
	[[file]|results].out:	This file contains a summary with pdb pairs superposed
		and RMSD values obtained by click
"""


def add_superposition_pairs(inputs, superposition_pairs):
    '''Read each line of superposition list file given by the user and creates
	two PDB objects spliting lines by comma ",". This objects will be appended
	to superposition_pair list var'''

    folder_path = ''

    for i in open(inputs.superpositions_list):
        if i[:5] == "PATH:":
            folder_path = i[5:].strip()
            if folder_path[-1] != '/':
                folder_path += '/'
        aux = i.split(',')

        if len(aux) == 2:
            superposition_pairs.append((PDB(folder_path + aux[0].strip()), \
                                        PDB(folder_path + aux[1].strip())))


def leer(archivo):
    f = open(archivo)
    num = 0
    rmsd = None
    for line in f:
        if num == 0:
            # ~ print line
            if line[28:].strip() == '0':
                # if this is iqual to 0 that means matched atoms was found :(
                rmsd = "----"
                break
            else:
                num += 1
                continue
        elif num == 1:
            try:
                rmsd = line[6:].strip()
                f.close()
                return rmsd
            except:
                f.close()
                return '----'
        num += 1
    f.close()


def click(superposition_list, flags, d_table, d_pdbResnumTable, d_table_inverse):
    '''This function takes a superposition_list of PDB objects and superpose
	them with click, then returns a results list with each superposition pair
	pdbs, RMSD and SO reported by click'''

    results = []
    # ~ fbasura = open("basura.basura","w")

    for mobile, target in superposition_list:
        if check(mobile["pdb"]) and check(target["pdb"]):
            # running click
            if '-s' in flags:
                # with this flag click doesn't creates rotated pdb files
                proc = subprocess.Popen(["basura", "%sMod.pdb" % mobile["pdb"][:-4], \
                                         "%sMod.pdb" % target["pdb"][:-4], \
                                         "-s 0"], \
                                        stdout=subprocess.PIPE,
                                        shell=False,
                                        executable='./click',
                                        close_fds=True,
                                        ).communicate()[0]
            # ~ proc = subprocess.call(["basura","%sMod.pdb"%mobile["pdb"][:-4],
            # ~ "%sMod.pdb"%target["pdb"][:-4],
            # ~ "-s 0"],
            # ~ shell = False,
            # ~ executable='./click',
            # ~ close_fds=True,
            # ~ )

            # ~ (mobile["pdb"][:-4],target["pdb"][:-4])],stdout=subprocess.PIPE ,shell = True,close_fds=True)
            else:
                proc = subprocess.Popen(["basura", "%sMod.pdb" % mobile["pdb"][:-4], \
                                         "%sMod.pdb" % target["pdb"][:-4], ],
                                        stdout=subprocess.PIPE,
                                        shell=False,
                                        executable='./click',
                                        close_fds=True,
                                        ).communicate()[0]

            lines = proc.split("\n")
            # ~ print mobile["name"],target["name"]
            # ~ print lines[0]

            if lines[0][28:].strip() == '0':
                # if this is iqual to 0 that means matched atoms was found :(
                rmsd = " NA "
            else:
                try:
                    rmsd = lines[1][6:].strip()
                except:
                    rmsd = " NA "
            results.append((mobile["name"], target["name"], rmsd))
            sys.stdout.flush()

            # ~ rmsd = leer("%s%sMod-%sMod.pdb.1.clique"%(mobile["path"], mobile["name"],target["name"]))
            # ~ results.append((mobile["name"],target["name"],rmsd))

            if not '-s' in flags:
                # merging click PDB files outputs
                mergeClickResults(flags, mobile, target)

                # modify the resulting PDB file to IUPAC normal nomenclature
                merged_pdb_full = "%s%s_%s_Merged.pdb" % (mobile["path"], mobile["name"], target["name"])
                merged_pdb_name = "%s_%s_Merged" % (mobile["name"], target["name"])
                try:
                    modifyPDBAtoms(flags, None, d_table_inverse, d_pdbResnumTable, merged_pdb_full, merged_pdb_name,
                                   clickOutput=True)
                except:
                    pass


        else:
            pdb1 = mobile["name"]
            pdb2 = target["name"]
            if not check(mobile["pdb"]):
                pdb1 = "(%s)" % mobile["name"]
            if not check(target["pdb"]):
                pdb2 = "(%s)" % target["name"]

            results.append((pdb1, pdb2, " NA "))

        # deleting all temp files created
        deleteModifiedFiles(flags, mobile, target)
    # eliminamos los archivos .clique creados
    # ~ if not '-m' in flags:
    # ~ proc = subprocess.Popen(["rm -f %s/%s-%s.pdb.1.clique"%\
    # ~ (mobile["path"], mobile["name"],target["name"])],stdout=subprocess.PIPE ,shell = True,close_fds=True)
    # ~ l2 = proc.communicate()
    # ~ fbasura.close()
    return results


def divide_lists(superposition_list, parts):
    '''returns "parts" lists of superposition pairs for each CPU'''
    rango = len(superposition_list) // parts + 1
    cuts = [i * rango for i in xrange(0, parts)]

    lists = []
    k = 1
    while k < len(cuts):
        lists.append(superposition_list[cuts[k - 1]:cuts[k]])
        k += 1
    lists.append(superposition_list[cuts[-1]:])

    return lists


def write_results(results, f_out):
    localtime = time.asctime(time.localtime(time.time()))
    f_out.write("End time :\t%s\n" % localtime)
    f_out.write("\n%s%s%s\n" % ("PDB1".ljust(30), "PDB2".ljust(30), "RMSD".rjust(10)))
    f_out.write("-" * 70 + "\n")
    for i in results:
        pdb1 = i[0]
        pdb2 = i[1]
        rmsd = i[2]
        f_out.write("%s%s%s\n" % (pdb1.ljust(30), pdb2.ljust(30), rmsd.rjust(10)))


######################
## PROGRAM COURSE ##
######################

if __name__ == '__main__':

    flags = sys.argv
    # Cheking inputs vars
    inputs, superposition_pairs = prepareInputs(flags)

    # saving tables info
    d_table, d_table_inverse, d_pdbResnumTable = tables(inputs)

    # opening file for write outputs
    if '-o' in flags:
        f_out = open(inputs.output, "w")
    else:
        f_out = open("results.output", "w")

    localtime = time.asctime(time.localtime(time.time()))
    f_out.write("Init time :\t%s\n" % localtime)

    # Changing Parameters.inp file according the table
    if not '-i' in flags:
        modifyClickParameters()
    # if -i flag is used, we change parameters.inp for the new parameters file
    else:
        replaceParametersInp(inputs.parameters)

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
