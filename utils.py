import sys

def helpInfo():
    print """
DEPENDENCIES

	This python script needs the next files in the same script folder:
		*Click binary
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

def errorSentence():
    # ~ print "\nTry one of this two options:\n1)python modclick pdb1.pdb pdb2.pdb -t tableFile"
    # ~ print "2)python modclick pdb1.pdb pdb2.pdb -t tableFile -i ParametersFile\n"
    helpInfo()
    sys.exit(1)