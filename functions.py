import sys
import subprocess
import os
import time
from utils import *
from functions import *
import cPickle
import Inputs
import PDB


###############
## FUNCTIONS ##
###############

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
        pdb = pdb_obj.pdb
    if pdbname:
        name_pdb = pdbname
    else:
        name_pdb = pdb_obj.name

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


def modifyPDBAtoms2(inputs, pdb_obj, fullpdb=None, pdbname=None, clickOutput=False):
    '''Take a pdb file as input and returns a new PDB file with atomtypes
	modifications. The new file has the "Mod" word before .pdb extension.
	If clickOutput flag is True, the pdb file will be modified with normal
	IUPAC names for a standard PDB file'''
    d_keywordRecord = {"ATOM": None, "HETATM": None, "ANISOU": None, "SIGUIJ": None}

    if fullpdb:
        pdb = fullpdb
    else:
        pdb = pdb_obj.pdb
    if pdbname:
        name_pdb = pdbname
    else:
        name_pdb = pdb_obj.name

    l_newPDB = []
    for line in open(pdb):
        if line[:6].strip() in d_keywordRecord:
            resname = line[17:20].strip()
            atomname = line[12:16].strip()
            chain = line[21]
            resnum = line[22:26].strip()

            # first we can modify a line just if the table has the right resname and atomname
            if inputs.d_table.has_key(resname):
                if inputs.d_table[resname].has_key(atomname):
                    if inputs.f:
                        # if -f flag is given, we use only the line with residue specified in d_pdbResnumTable
                        if name_pdb in inputs.d_pdbResnumTable and \
                           chain in inputs.d_pdbResnumTable[name_pdb] and \
                           resnum in inputs.d_pdbResnumTable[name_pdb][chain]:

                            l_newPDB.append(line[:12] + formatAtomName(inputs.d_table[resname][atomname]) + line[16:])
                        else:
                            # ~ l_newPDB.append(line[:12] + formatAtomName(d_table[resname][atomname]) + line[16:])
                            l_newPDB.append(line)
                    else:
                        l_newPDB.append(line[:12] + formatAtomName(inputs.d_table[resname][atomname]) + line[16:])
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



def deleteModifiedFiles(inputs, mobile, target, deleteModfiles=False):
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
            os.unlink("%sMod.pdb" % (mobile.pdb[:-4]))
            os.unlink("%sMod.pdb" % (target.pdb[:-4]))
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
            os.unlink("%sMod-%sMod.1.pdb" % (mobile.pdb[:-4], target.name))
            os.unlink("%sMod-%sMod.1.pdb" % (target.pdb[:-4], mobile.name))
        except:
            pass

        # if -m flag is not given, then try to delete generated .clique files
        if not inputs.m:
            # ~ proc = subprocess.call(["rm -f %s%sMod-%sMod.pdb.1.clique"%(mobile["path"],
            # ~ mobile["name"],
            # ~ target["name"])],
            # ~ shell = True)
            try:
                os.unlink("%s%sMod-%sMod.pdb.1.clique" % (mobile.path,
                                                          mobile.name,
                                                          target.name))
                os.unlink("%s%s_%s_Merged.pdb" % (mobile.path, mobile.name, target.name))
            except:
                pass


def mergeRotatedMobileWithTarget(inputs, mobile, target):
    '''Merge rotated mobile click result with target in a single file. The files are separated with
	MODEL record'''

    try:
        f_mobile = open("%s%sMod-%sMod.1.pdb" % (mobile.path, mobile.name, target.name))

        # ~ f3 = open("%s%sMod-%sMod.1.pdb"%(target["path"], target["name"], mobile["name"]))
        # target doesn't change its coordinates
        # f_target = open("%s%sMod.pdb" % (target.path, target.name))

    except:
        deleteModifiedFiles(inputs, mobile, target)
        # ~ sys.exit(1)
        return 0

    # Open merged file to write pdb lines
    f = open("%s%s_%s_Merged.pdb" % (mobile.path, mobile.name, target.name), "w")
    f.write("%s\n" % ("MODEL        1".ljust(80)))


    # print "largo modificado"
    # print len(f_mobile.readlines())
    # print "largo original"
    # print len(mobile.original_lines)

    wa = f_mobile.readlines()
    # print "last line"
    # print wa[-1]
    for i,j in enumerate(wa):
        # if i== 0:
            # print j
            # print mobile.original_lines[i]

        # if len(wa) == i+2:
            # print j
            # print mobile.original_lines[i]

        # print j
        if len(j.strip()) > 0:
            coordinate_section = j[30:54]
            modified_line = mobile.original_lines[i][:30] + coordinate_section + mobile.original_lines[i][54:]
            # print modified_line
            f.write(modified_line)

    f.write("%s\n" % "ENDMDL".ljust(80))
    f.write("%s\n" % ("MODEL        2".ljust(80)))

    # for i in f_target.readlines():
    for i in target.original_lines:
        f.write(i)
    f.write("ENDMDL".ljust(80))

    f.close()
    f_mobile.close()


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
            superposition_pairs.append((PDB.PDB(folder_path + aux[0].strip(), inputs),
                                        PDB.PDB(folder_path + aux[1].strip(), inputs)))


# def leer(archivo):
#     f = open(archivo)
#     num = 0
#     rmsd = None
#     for line in f:
#         if num == 0:
#             # ~ print line
#             if line[28:].strip() == '0':
#                 # if this is iqual to 0 that means matched atoms was found :(
#                 rmsd = "----"
#                 break
#             else:
#                 num += 1
#                 continue
#         elif num == 1:
#             try:
#                 rmsd = line[6:].strip()
#                 f.close()
#                 return rmsd
#             except:
#                 f.close()
#                 return '----'
#         num += 1
#     f.close()


def click(superposition_list, inputs):
    '''This function takes a superposition_list of PDB objects and superpose
	them with click, then returns a results list with each superposition pair
	pdbs, RMSD and SO reported by click'''


    # clickExecutable = './click'
    #clickExecutable = '/home/norodrig/Programs/bin/Click/click'
    clickExecutable = '/home/norodrig/Programs/bin/Click_multiple/clickmul'

    results = []
    # ~ fbasura = open("basura.basura","w")

    for mobile, target in superposition_list:

        itsOk1 = check(mobile.pdb)
        itsOk2 = check(target.pdb)

        # if check(mobile.pdb) and check(target.pdb):
        if itsOk1 and itsOk2:

            # running click
            # if -s flag was given
            if inputs.s:
                # print "%sMod.pdb" % mobile.pdb[:-4]
                # print "%sMod.pdb" % target.pdb[:-4]
                commands1 = ["basura", "%sMod.pdb" % mobile.pdb[:-4],
                            "%sMod.pdb" % target.pdb[:-4],
                            "-s 0"
                            ]

                commands2 = ["basura", "%sMod.pdb" % target.pdb[:-4],
                            "%sMod.pdb" % mobile.pdb[:-4],
                            "-s 0"
                            ]

            else:

                commands1 = ["basura", "%sMod.pdb" % mobile.pdb[:-4],
                            "%sMod.pdb" % target.pdb[:-4]
                            ]

                commands2 = ["basura", "%sMod.pdb" % target.pdb[:-4],
                            "%sMod.pdb" % mobile.pdb[:-4]
                            ]


            # Run click process
            proc1 = subprocess.Popen(commands1,
                                    stdout=subprocess.PIPE,
                                    shell=False,
                                    executable=clickExecutable,
                                    close_fds=True,
                                    ).communicate()[0]

            proc2 = subprocess.Popen(commands2,
                        stdout=subprocess.PIPE,
                        shell=False,
                        executable=clickExecutable,
                        close_fds=True,
                        ).communicate()[0]

            #print proc1, proc2
            # Now extract RMSD value for the records
            lines1 = proc1.split("\n")
            lines2 = proc2.split("\n")
            print lines1
            print lines2

            for l, ll in zip(lines1, lines2):
                if "The number of matched residues" in l:
                    match1 = int(l.split(" ")[-1])
                if "The number of matched residues" in ll:
                    match2 = int(ll.split(" ")[-1])
                if "RMSD" in l:
                    rmsd1 = float(l.split(" ")[-1])
                if "RMSD" in ll:
                    rmsd2 = float(ll.split(" ")[-1])
                if "Structure Overlap" in l:
                    strover1 = float(l.split(" ")[-1])
                if "Structure Overlap" in ll:
                    strover2 = float(ll.split(" ")[-1])


            if rmsd1 < rmsd2 and match2 <= match1:
                rmsd = rmsd1
                print mobile.name + "/" + target.name
                print "RMSD: " + str(rmsd) + " (Matches: " + str(match1) + " , Structural Overlap: " + str(strover1) +")"
                results.append((mobile.name, target.name, rmsd1, match1, strover1))
            else:
                rmsd = rmsd2
                print target.name + "/" + mobile.name
                print "RMSD: " + str(rmsd) + " (Matches: " + str(match2) + " , Structural Overlap: " + str(strover2) +")"
                results.append((mobile.name, target.name, rmsd2, match2, strover2))

            sys.stdout.flush()


            # Remember that -s flag means runs click with -s flag. With that click doesn't create transformed pdbs
            # and cannot generate Merged pdb
            if not inputs.s:
                # print "%sMod.pdb" % mobile.pdb[:-4]
                # merging click PDB files outputs
                mergeRotatedMobileWithTarget(inputs, mobile, target)


        else:
            if itsOk1:
                name1 = mobile.name
            else:
                name1 = "(%s)"%mobile.name

            if itsOk2:
                name2 = target.name
            else:
                name2 = "(%s)"%target.name
            results.append((name1, name2, '----'))
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
    f_out.write("\n%s%s%s%s\n" % ("PDB1".ljust(40), "PDB2".ljust(40), "RMSD".rjust(10), "  MatchRes/StrOverlap".rjust(10)))
    f_out.write("-" * 70 + "\n")
    for i in results:
        pdb1 = i[0]
        pdb2 = i[1]
        rmsd = str(i[2])
        match = str(i[3])
        strover = str(i[4])
        f_out.write(pdb1+"\t"+pdb2+"\t"+rmsd+"\t"+match+"\t"+strover+"\n")