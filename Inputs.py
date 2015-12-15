from utils import *
import sys
import PDB
import functions
import subprocess
import os


class Inputs:
    '''Inputs class creates a object with command line entry inputs'''

    def __init__(self, args):
        self.parameters = None
        self.d_table = {}  # dictionary with {atomName:newAtomType}
        self.tablefile = None  # file given by the user that contain list of new atomtypes
        self.d_pdbResnumTable = {}  # dictionary with {PDB:{CHAIN:{resnum:None,},},}
        self.pdb_resnum_table_file = None
        self.superpositions_list = None  # superposition list file if given
        self.numCPUs = None
        self.output = None
        self.original_PDB_lines = []  # with this I can return pdb to its IUPAC atom names convention
        self.pdb1 = None # pdb file name given by the user
        self.pdb2 = None # pdb file name given by the user
        # self.superposition_pairs = []  # list of all superposition pairs

        # INPUTS FLAGS

        # -o [file] modclick writes summary [file] with pdb1,pdb2,RMSD records. If not
        # [file] name is given or this flag is not used, results.out file will be
        # written
        self.o = False

        # -i [file] modclick uses the [file] parameters to run click superpositions
        self.i = False

        # -f [file] modclick uses the [file] table with PDB:chain:resnumber fields that
		# shows only the residues wich you want to change atom types
        self.f = False

        # -p [nCPUs] modclick try to run superpositions with parallel python module using
		# a given number of cores. If no number is specified, or 0 number is given,
		# modclick uses the maximun machine available cores
        self.p = False

        # -s	modclick runs click with -s 0 flag. With this flag click doesn't create
		# transformation pdb files, therefore Merged pdb file is not created
        self.s = False

        # -m	modclick mantains *.clique file. Merged pdb file is mantains too except if
		# -s flag is chosen
        self.m = False

        self.nom_atomos = {
            'A': ['N1', 'C2', 'N3', 'C4', 'C5', 'C6', 'N7', 'C8', 'N9', 'N6'],
            'G': ['N1', 'C2', 'N3', 'C4', 'C5', 'C6', 'N7', 'C8', 'N9', 'O6', 'N2'],
            'C': ['N1', 'C2', 'N3', 'C4', 'C5', 'C6', 'O2', 'N4'],
            'T': ['N1', 'C2', 'N3', 'C4', 'C5', 'C6', 'O2', 'O4', 'C7'],
            # notar que el C7 es el metilo en el C5 \
            'U': ['N1', 'C2', 'N3', 'C4', 'C5', 'C6', 'O2', 'O4']
        }

        self.__saveAndCheckInputs(args)
        self.__extractInfoTables()
        self.__checkParametersInp()


    def __saveAndCheckInputs(self, args):
        """With this function we iterate for all args and saves the inputs"""

        # creating inputs object
        # inputs = Inputs()
        # superposition_pairs = []  # list of all superposition pairs
        par = []  # var use for save pdb inputs


        # first iterate fdor all argv and recopilate flag information
        for j, i in enumerate(args):
            # save pdf file if exist
            if i.strip()[-4:] == '.pdb':
                if self.pdb1 == None:
                    self.pdb1 = i.strip()
                # par.append(i.strip())
                else:
                    self.pdb2 = i.strip()
            # checking parameters file
            if i.strip() == '-i':
                try:
                    self.i = True
                    self.parameters = args[j + 1]
                except:
                    print "No Parameters file indicated"
                    errorSentence()

            # checking table file
            if i.strip() == '-t':
                try:
                    self.tablefile = args[j + 1]
                except:
                    print "No table file indicated"
                    errorSentence()

            # checking pdb_resnum_table
            if i.strip() == '-f':
                try:
                    self.f = True
                    self.pdb_resnum_table = args[j + 1]
                except:
                    print "No PDB_resnum_table file indicated"
                    errorSentence()

            # checking superposition list
            if i.strip() == '-l':
                try:
                    self.superpositions_list = args[j + 1]
                except:
                    print "No superpositions list file indicated"
                    errorSentence()

            # checking parallel python nCPUs
            if i.strip() == '-p':
                self.p = True
                try:
                    if args[j + 1].isdigit():
                        self.numCPUs = int(args[j + 1])
                except:
                    pass

            # checking if output file name is given
            if i.strip() == '-o':
                try:
                    self.o = True
                    self.output = args[j + 1].strip()
                except:
                    pass

            # check -s flag. If true, click runs with -s flag
            if i.strip() == '-s':
                self.s = True

            # check -m flag. If true, mantains *.clique file
            if i.strip() == '-m':
                self.m = True

        # cheking if the files exists
        if not len(par) in (0, 2):
            print "two pdf files are necessary as minimum or a list of superposition pairs."
            errorSentence()

        # if par and par[0] and par[1]:
        if self.pdb1 and self.pdb2:
            # if user gives two pdb files
            if not functions.check(self.pdb1):
                print "%s does'n exist, please try a real pdb file." % self.pdb1
                errorSentence()

            if not functions.check(self.pdb2):
                print "%s does'n exist, please try a real pdb file." % self.pdb2
                errorSentence()

        else:
            # if no list of superposition is given
            if not '-l' in args:
                print "Click needs two pdb files."
                errorSentence()

        # if check
        if not functions.check(self.parameters):
            print "%s does'n exist, please try a real parameters file." % self.parameters
            errorSentence()

        if not self.tablefile:
            print "modclick needs a table file."
            errorSentence()
        else:
            if not functions.check(self.tablefile):
                print "%s does'n exist, please try a real table file." % self.tablefile
                errorSentence()

        if not functions.check(self.pdb_resnum_table_file):
            print "%s does'n exist, please try a real PDB_resnum_table file." % self.pdb_resnum_table_file
            errorSentence()

        # checking superpsition list file
        if not functions.check(self.superpositions_list):
            print "%s does'n exist, please try a real superposition list." % self.superpositions_list
            errorSentence()

        # checking Parameters.inp file
        if not functions.check("Parameters.inp"):
            print "Parameters.inp does'n exist."
            sys.exit(1)

    def __extractInfoTables(self):
        '''Open the "table" file with new atomtypes specifications and save
        the info into a dicctionary "d_table". A inverse dictionary "d_table_inverse"
        is also created to restore PDB files to IUPAC nomenclature.'''

        # d_table_inverse = {}  # dictionary inverse to d_table
        for i, line in enumerate(open(self.tablefile)):
            if i > 0 and len(line.strip()) > 0:
                aux = line.split(":")
                if len(aux) == 3:
                    if not self.d_table.has_key(aux[0].strip()):
                        self.d_table[aux[0].strip()] = {aux[1].strip(): aux[2].strip(), }
                        # d_table_inverse[aux[0].strip()] = {aux[2].strip(): aux[1].strip(), }
                    else:
                        self.d_table[aux[0].strip()][aux[1].strip()] = aux[2].strip()
                        # d_table_inverse[aux[0].strip()][aux[2].strip()] = aux[1].strip()

        # now if -f flag is choosen
        if self.pdb_resnum_table_file:
            f = open(self.pdb_resnum_table_file)
            l = f.readlines()
            f.close()
            for i in l[1:]:
                aux = i.split()
                pdb = aux[0].strip()
                chain = aux[1].strip()
                number = aux[2].strip()
                if not self.d_pdbResnumTable.has_key(pdb):
                    self.d_pdbResnumTable[aux[0]] = {}
                if not self.d_pdbResnumTable[aux[0]].has_key(chain):
                    self.d_pdbResnumTable[pdb][chain] = {}
                self.d_pdbResnumTable[pdb][chain][number] = None
        else:
            self.d_pdbResnumTable = None

            # return d_table, d_table_inverse, d_pdbResnumTable

    def __modifyParametersInp(self):
        '''Modify the Parameters.inp file according the table information'''

        # Creating new typeAtom line
        newTypeAtomLine = "typeAtom=\""
        typeatoms = []
        for i in self.d_table:
            typeatoms.extend(self.d_table[i].values())
        daux = {i: None for i in typeatoms}
        typeatoms = daux.keys()
        typeatoms.sort()
        for i, j in enumerate(typeatoms):
            if i == len(typeatoms) - 1:
                newTypeAtomLine += functions.formatAtomName(j)
            else:
                newTypeAtomLine += "%s ," % functions.formatAtomName(j)
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


    def replaceParametersInp(self, restore=False):
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
            for i in open(self.parameters):
                f.write(i)
            f.close()
        else:
            f = open("Parameters.inp", "w")
            for i in open("Parameters.inp.backup"):
                f.write(i)
            f.close()

            #removes the backpup
            try:
                os.remove('Parameters.inp.backup')
            except OSError:
                pass
            # proc = subprocess.Popen(["rm -f Parameters.inp.backup"], stdout=subprocess.PIPE, shell=True)
            # l2 = proc.communicate()


    def __checkParametersInp(self):
        # if -i flag is used, we change parameters.inp for the new parameters file
        if self.i:
            self.replaceParametersInp(self.parameters)
        else:
            # Changing Parameters.inp file according the table
            self.__modifyParametersInp()

